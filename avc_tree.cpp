/*/
Copyright © 2011 DigiPen Institute of Technology

Author: Joel Barba
Date: 11-27
Purpose: Implements an AVC tree.

/*/

#include "avc_tree.hpp"
#include "random.hpp"
#include "axis_aligned_plane.hpp"
#include "box_extents.hpp"

#include <ppl.h>

namespace
{
	namespace sides
	{
		enum type
		{
			top,
			bottom,
			left,
			right
		};
	}

	std::unique_ptr<avc_node> build_tree(Box3D pVolume, avc_node::segment_vector pSegments,
		avc_node::potentially_visible_set pPVS, const avc_node::potentially_visible_set& pMainPVS, 
		float pCost, unsigned pDepth = 0);

	bool is_leaf_node(unsigned pDepth)
	{
		return pDepth == 6;
	}

	Vector3D create_point(std::vector<unsigned>& pSides)
	{
		Vector3D point(1, 1, 0.01f);
		auto indexA = 0, indexB = 1;

		unsigned index = rand() % pSides.size();
		auto side = pSides[index];
		pSides.erase(pSides.begin() + index);

		if(side == sides::left || side == sides::right)
			std::swap(indexA, indexB);

		if(side == sides::left || side == sides::bottom)
			point[indexB] = -1;

		point[indexA] = homogenous_random() * 2 - 1;

		return point;
	}

	view_segment create_max_view_segment()
	{
		std::vector<unsigned> sideVector;
		sideVector.push_back(sides::top);
		sideVector.push_back(sides::bottom);
		sideVector.push_back(sides::left);
		sideVector.push_back(sides::right);

		Segment3D segment(create_point(sideVector), create_point(sideVector));

		return view_segment(std::move(segment));
	}

	float avc_heuristic(const box_extents& pVolume, unsigned pPvsSize)
	{
		return volume(pVolume) * pPvsSize;
	}

	void choose_split(axis_aligned_plane& pPlane, const Box3D& pVolume)
	{
		pPlane.axis = rand() % 2;
		pPlane.intersection = pVolume.center[pPlane.axis];
	}

	std::pair<std::unique_ptr<avc_node>, std::unique_ptr<avc_node>>
		dumb_split(Box3D pVolume, avc_node::segment_vector pSegments, avc_node::potentially_visible_set pPVS,
		const avc_node::potentially_visible_set& pMainPVS, float pCost, unsigned pDepth)
	{
		// find the potential split plane.
		auto axis = pDepth % 2;
		axis_aligned_plane plane(axis, pVolume.center[axis]);

		// split the parent volume by the plane.
		box_extents leftVolume, rightVolume;
		split(leftVolume, rightVolume, pVolume, plane);

		// split the parent segments by the plane.
		avc_node::segment_vector leftSegments, rightSegments;
		split(leftSegments, rightSegments, pSegments, plane);

		// calculate the sizes of the left and right PVS.
		// calculate the sizes of the left and right PVS.
		avc_node::potentially_visible_set leftPVS, rightPVS;
		Concurrency::parallel_invoke([&]{ leftPVS = split(pMainPVS, leftSegments); },
			[&]{ rightPVS = split(pMainPVS, rightSegments); });

		// calculate the cost.
		auto leftCost = avc_heuristic(leftVolume, leftPVS.size());
		auto rightCost = avc_heuristic(rightVolume, rightPVS.size());

		// build child trees.
		auto left = build_tree(to_box(leftVolume), std::move(leftSegments), std::move(leftPVS),
			pMainPVS, leftCost, pDepth + 1);
		auto right = build_tree(to_box(rightVolume), std::move(rightSegments), std::move(rightPVS),
			pMainPVS, rightCost, pDepth + 1);

		return std::make_pair(std::move(left), std::move(right));
	}

	std::pair<std::unique_ptr<avc_node>, std::unique_ptr<avc_node>>
		split(Box3D pVolume, avc_node::segment_vector pSegments, avc_node::potentially_visible_set pPVS,
				const avc_node::potentially_visible_set& pMainPVS, float pCost, unsigned pDepth)
	{
		box_extents bestLeftVolume, bestRightVolume;
		avc_node::segment_vector bestLeftSegments, bestRightSegments;
		avc_node::potentially_visible_set bestLeftPVS, bestRightPVS;
		float bestLeftCost = 0, bestRightCost = 0, bestCost = std::numeric_limits<float>::max();
		axis_aligned_plane bestPlane;
		
		unsigned i = 0;
		for (auto viewSegment = pSegments.begin(); viewSegment != pSegments.end() && i < 12; ++viewSegment, ++i)
		{
			for(unsigned segmentPoint = 0; segmentPoint < 2; ++segmentPoint)
			{
				for (unsigned axis = 0; axis < 2; ++axis)
				{
					// find the potential split plane.
					axis_aligned_plane plane(axis, viewSegment->segment[segmentPoint][axis]);

					// split the parent volume by the plane.
					box_extents leftVolume, rightVolume;
					split(leftVolume, rightVolume, pVolume, plane);

					// split the parent segments by the plane.
					avc_node::segment_vector leftSegments, rightSegments;
					split(leftSegments, rightSegments, pSegments, plane);

					// calculate the sizes of the left and right PVS.
					avc_node::potentially_visible_set leftPVS, rightPVS;
					Concurrency::parallel_invoke([&]{ leftPVS = split(pMainPVS, leftSegments); },
													[&]{ rightPVS = split(pMainPVS, rightSegments); });

					// calculate the cost.
					auto leftCost = avc_heuristic(leftVolume, leftPVS.size());
					auto rightCost = avc_heuristic(rightVolume, rightPVS.size());

					auto cost = leftCost + rightCost - pCost;

					//set the best split plane.
					if(bestCost > cost)
					{
						bestLeftVolume = leftVolume;
						bestRightVolume = rightVolume;
						bestLeftSegments = leftSegments;
						bestRightSegments = rightSegments;
						bestLeftPVS = leftPVS;
						bestRightPVS = rightPVS;
						bestLeftCost = leftCost;
						bestRightCost = rightCost;
						bestCost = cost;
						bestPlane = plane;
					}
				}
			}
		}

		// build child trees.
		auto left = build_tree(to_box(bestLeftVolume), std::move(bestLeftSegments), std::move(bestLeftPVS),
								pMainPVS, bestLeftCost, pDepth + 1);
		auto right = build_tree(to_box(bestRightVolume), std::move(bestRightSegments), std::move(bestRightPVS),
								pMainPVS, bestRightCost, pDepth + 1);

		return std::make_pair(std::move(left), std::move(right));
	}

	std::unique_ptr<avc_node> build_tree(Box3D pVolume, avc_node::segment_vector pSegments,
		avc_node::potentially_visible_set pPVS, const avc_node::potentially_visible_set& pMainPVS, 
		float pCost, unsigned pDepth)
	{
		// terminate with leaf node.
		if(is_leaf_node(pDepth))
			return std::unique_ptr<avc_node>(new avc_node(std::move(pVolume),
															std::move(pSegments), std::move(pPVS)));

		// find the best splitting plane.
		auto children = split(std::move(pVolume), std::move(pSegments), std::move(pPVS),
			pMainPVS, pCost, pDepth);

		// construct internal node with children.
		return std::unique_ptr<avc_node>(new avc_internal_node(std::move(children.first),
																std::move(children.second)));
	}

	const avc_node::potentially_visible_set* traverse_tree(const std::unique_ptr<avc_node>& pNode,
															Vector3D pPosition)
	{
		const avc_node::potentially_visible_set* pvs = nullptr;

		if(pNode->volume.contains(pPosition))
		{
			pvs = &pNode->pvs;
		}
		else
		{
			if (pNode->left)
				pvs = traverse_tree(pNode->left, pPosition);
			if (!pvs && pNode->right)
				pvs = traverse_tree(pNode->right, pPosition);
		}

		return pvs;
	}
}

avc_tree::avc_tree()
{
	srand(100);
}

void avc_tree::build(const Box3D& pVolume, unsigned pRayCount, const std::vector<Box3D>& pObjects)
{
	//*/
	// generate the maximum view segments.
	avc_node::segment_vector segments;
	segments.reserve(pRayCount);
	for(unsigned i = 0; i < pRayCount; ++i)
		segments.push_back(create_max_view_segment());
	/*/
	// generate the maximum view segments.
	avc_node::segment_vector segments;
	segments.push_back(Segment3D(Vector3D(1, 1, 0.1f), Vector3D(-1, -1, 0.1f)));
	segments.push_back(Segment3D(Vector3D(-1, 0.25f, 0.1f), Vector3D(.3f, 1, 0.1f)));
	segments.push_back(Segment3D(Vector3D(-1, -.25f, 0.1f), Vector3D(1, .5f, 0.1f)));
	//*/

	// generate the potentially visible set.
	avc_node::potentially_visible_set pvsMain;
	for(auto object = pObjects.begin(); object != pObjects.end(); ++object)
		pvsMain.insert(&*object);

	// find the PVS of the root node.
	auto pvs = split(pvsMain, segments);

	//calculate the cost of the node.
	auto cost = avc_heuristic(pVolume, pvs.size());

	mRoot = build_tree(pVolume, std::move(segments), std::move(pvs), pvsMain, cost);
}

const avc_node::potentially_visible_set* avc_tree::traverse(Vector3D& pPosition) const
{
	if (!mRoot)
		return nullptr;
		
	return traverse_tree(mRoot, pPosition);
}
