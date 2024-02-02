# Octree/Quadtree/N-dimensional linear tree
[![MSBuild and Unittests](https://github.com/attcs/Octree/actions/workflows/msbuild.yml/badge.svg)](https://github.com/attcs/Octree/actions/workflows/msbuild.yml)
<br>
Lightweight, parallelizable C++ implementation of an Octree/Quadtree/N-d orthotree using Morton Z curve-based location code ordering.<br>
<br>
What is the Octree and what is good for? https://en.wikipedia.org/wiki/Octree

## Features
* Adaptable to any existing geometric system
* Arbitrary number of dimensions for other scientific usages
* Support of `std::execution` policies (so it is parallelizable)
* Edit functions to Insert/Update/Erase entities
* Wide range of search functions
  * Range search
  * Pick search
  * K - Nearest neighbor search
  * Ray-traced search
* Collision detection
* Nodes can be accessed in O(1) time
* Search is accelerated by Morton Z curve based location code
* Both the non-owning `Core` and the `Container` wrapper is provided

## Limitations
* Maximum number of dimensions is 63.
* Maximum depth of octree solutions is 10.
* Abstract classes cannot be used for `vector_type` and `box_type`

## Requirements
* Language standard: C++20 or above

## Time complexity
* Creation: O(n)
* Range search: O(log{2^N}(n)) (where N is the number of the dimension)
* Access any node: O(1)

## Usage
* Use `AdaptorBasicsConcept` or `AdaptorConcept` to adapt the actual geometric system. It is not a necessary step, basic point/vector and bounding box objects are available.
* Call the static member function `Create()` for a contiguous container (any `std::span` compatible) of Points or Bounding boxes to build the tree. It supports `std::execution` policies (e.g.: `std::execution::parallel_unsequenced_policy`) which can be effectively used to parallelize the creation process. (Template argument of the `Create()` functions)
* Call `PickSearch()` / `RangeSearch()` member functions to collect the wanted id-s
* Call `Core` edit functions `Insert()`, `Update()`, `UpdateIndexes()`, `Erase()` if the some of the underlying geometrical elements were changed or reordered
* Call `Container` edit functions `Add()`, `Update()`, `Erase()` if one of the underlying geometrical element was changed 
* Call `CollisionDetection()` member function for bounding box overlap examination.
* Call `VisitNodes()` to traverse the tree from up to down (breadth-first search) with user-defined `selector()` and `procedure()`.
* Call `GetNearestNeighbors()` for kNN search in point based tree. https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm
* Call `RayIntersectedFirst()` or `RayIntersectedAll()` to get intersected bounding boxes in order by a ray.


## Notes
* Header only implementation.
* Point and Bounding box-based solution is distinguished.
* Core types store only the entity ids, use Container types to store. Core types advantages: not copying and managing the entity information; disadvantages: this information may have to be provided again for the member function call.
* Container types have "C" postfix (e.g.: core `OctreeBox`'s container is `OctreeBoxC`).
* Bounding box-based solution stores item id in the parent node if it is not fit into any child node. Using `nSplitStrategyAdditionalDepth` template parameter, these boxes can be splitted then placed on the deeper level of the tree. The `nSplitStrategyAdditionalDepth` default is 2 and this split method is applied by default.
* Edit functions are available but not recommended to majorly build the tree.
* If less element is collected in a node than the max element then the child node won't be created.
* The underlying container is a hash-table (`std::unordered_map`) under 16D, which only stores the id-s and the bounding box of the child nodes.
* Original geometry data is not stored, so any search function needs them as an input.
* Unit tests are attached. (Microsoft Unit Testing Framework for C++)
* Tested compilers: MSVC 2019, Clang 12.0.0, GCC 11.3

## Major aliases in OrthoTree
```C++
  /// Default geometrical base elements

  using Point2D = OrthoTree::PointND<2>;
  using Point3D = OrthoTree::PointND<3>;
  using BoundingBox2D = OrthoTree::BoundingBoxND<2>;
  using BoundingBox3D = OrthoTree::BoundingBoxND<3>;


  /// Core types

  // Quadtree for points (2D)
  using QuadtreePoint = TreePointND<2>;

  // Quadtree for bounding boxes (2D)
  using QuadtreeBox = TreeBoxND<2>;

  // Octree for points (3D)
  using OctreePoint = TreePointND<3>;

  // Octree for bounding boxes (3D)
  using OctreeBox = TreeBoxND<3>;

  // Hexatree for points (4D)
  using HexatreePoint = TreePointND<4>;
  
  // ...
  using TreePoint16D = TreePointND<16>;


  /// Container types

  // Quadtree for points (2D)
  using QuadtreePointC = TreePointContainerND<2>;

  // Quadtree for bounding boxes (2D)
  using QuadtreeBoxC = TreeBoxContainerND<2>;

  // Octree for points (3D)
  using OctreePointC = TreePointContainerND<3>;

  // Octree for bounding boxes (3D)
  using OctreeBoxC = TreeBoxContainerND<3>;
```


## Basic examples

Usage of Container types
```C++
    #include "octree.h"
    using namespace OrthoTree;
    
    // Example #1: Octree for points
    {
      auto constexpr points = array{ Point3D{0,0,0}, Point3D{1,1,1}, Point3D{2,2,2} };
      auto const octree = OctreePointC(points, 3 /*max depth*/);

      auto const search_box = BoundingBox3D{ {0.5, 0.5, 0.5}, {2.5, 2.5, 2.5} };
      auto ids = octree.RangeSearch(search_box); // -> { 1, 2 }
      auto knn_ids = octree.GetNearestNeighbors(Point3D{ 1.1,1.1,1.1 }, 2 /*k*/); // -> { 1, 2 }
    }
    
    // Example #2: Quadtree for bounding boxes
    {
      auto boxes = vector
      {
        BoundingBox2D{ { 0.0, 0.0 }, { 1.0, 1.0 } },
        BoundingBox2D{ { 1.0, 1.0 }, { 2.0, 2.0 } },
        BoundingBox2D{ { 2.0, 2.0 }, { 3.0, 3.0 } },
        BoundingBox2D{ { 3.0, 3.0 }, { 4.0, 4.0 } },
        BoundingBox2D{ { 1.2, 1.2 }, { 2.8, 2.8 } }
      };

      auto quadtreebox = QuadtreeBoxC(boxes, 3
        , std::nullopt // user-provided bounding box for all
        , 2            // max element in a node 
        , false        // parallel calculation flag
      );

      // Collision detection
      auto ids_pairs_colliding = quadtreebox.CollisionDetection(); // { {1,4}, {2,4} }

      // Range search
      auto search_box = BoundingBox2D{ { 1.0, 1.0 }, { 3.1, 3.1 } };
      auto ids_inside = quadtreebox.RangeSearch(search_box); // -> { 1, 2, 4 }
      auto ids_overlaping = quadtreebox.RangeSearch<false /*overlap is enough*/>(search_box); 
        // -> { 1, 2, 3, 4 }

      // Picked boxes
      auto ptPick = Point2D{ 2.5, 2.5 };
      auto ids_picked = quadtreebox.PickSearch(ptPick); // -> { 2, 4 }
    }
    
    // Example #3: Parallel creation of octree for bounding boxes
    {
      auto boxes = vector{ BoundingBox3D{ { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 } } /* and more... */ };
      // Using ctor
      {
        auto octreebox = OctreeBoxC(boxes, 3, std::nullopt, OctreeBox::max_element_default
          , true // Set std::execution::parallel_unsequenced_policy
        );
      }
      // Using Create
      {
        auto octreebox = OctreeBoxC::Create<std::execution::parallel_unsequenced_policy>(boxes, 3);
      }
    }
```

Usage of Core types
```C++
    #include "octree.h"
    using namespace OrthoTree;
    
    // Example #1: Octree for points
    {
      auto constexpr points = array{ Point3D{0,0,0}, Point3D{1,1,1}, Point3D{2,2,2} };
      auto const octree = OctreePoint(points, 3 /*max depth*/);
       
      auto const search_box = BoundingBox3D{ {0.5, 0.5, 0.5}, {2.5, 2.5, 2.5} };
      auto ids = octree.RangeSearch(search_box, points); // -> { 1, 2 }
      auto knn_ids = octree.GetNearestNeighbors(Point3D{ 1.1,1.1,1.1 }, 2 /*k*/, points); // -> { 1, 2 }
    }
    
    // Example #2: Quadtree for bounding boxes
    {
      auto boxes = vector
      {
        BoundingBox2D{ { 0.0, 0.0 }, { 1.0, 1.0 } },
        BoundingBox2D{ { 1.0, 1.0 }, { 2.0, 2.0 } },
        BoundingBox2D{ { 2.0, 2.0 }, { 3.0, 3.0 } },
        BoundingBox2D{ { 3.0, 3.0 }, { 4.0, 4.0 } },
        BoundingBox2D{ { 1.2, 1.2 }, { 2.8, 2.8 } }
      };

      auto quadtreebox = QuadtreeBox(boxes, 3
        , std::nullopt // user-provided bounding box for all
        , 2            // max element in a node 
      );

      // Collision detection
      auto ids_pairs_colliding = quadtreebox.CollisionDetection(boxes); // { {1,4}, {2,4} }

      // Range search
      auto search_box = BoundingBox2D{ { 1.0, 1.0 }, { 3.1, 3.1 } };
      auto ids_inside = quadtreebox.RangeSearch(search_box, boxes); // -> { 1, 2, 4 }
      auto ids_overlaping = quadtreebox.RangeSearch<false/*overlap is enough*/>(search_box, boxes);
        // -> { 1, 2, 3, 4 }
      
      // Picked boxes
      auto ptPick = Point2D{ 2.5, 2.5 };
      auto ids_picked = quadtreebox.PickSearch(ptPick, boxes); // -> { 2, 4 }
    }
```
    
For more examples, see the unit tests.
<div align="center" width="100%"><img src="https://github.com/attcs/Octree/blob/master/docs/quadtree_example.PNG " align="center" height="300"></div>


## Adapting Octree/Quadtree to user-defined Point / Bounding box objects
```C++
  // User-defined geometrical objects

  struct Point2DCustom { float x; float y; };
  using BoundingBox2DCustom = std::array<Point2DCustom, 2>;


  // Adaptor

  struct AdaptorBasicsCustom
  {
    static inline float& point_comp(Point2DCustom& pt, OrthoTree::dim_type iDimension)
    {
      switch (iDimension)
      {
        case 0: return pt.x;
        case 1: return pt.y;
        default: assert(false); return pt.x;
      }
    }

    static constexpr float point_comp_c(Point2DCustom const& pt, OrthoTree::dim_type iDimension)
    {
      switch (iDimension)
      {
        case 0: return pt.x;
        case 1: return pt.y;
        default: assert(false); return pt.x;
      }
    }

    static inline Point2DCustom& box_min(BoundingBox2DCustom& box) { return box[0]; }
    static inline Point2DCustom& box_max(BoundingBox2DCustom& box) { return box[1]; }
    static constexpr Point2DCustom const& box_min_c(BoundingBox2DCustom const& box) { return box[0]; }
    static constexpr Point2DCustom const& box_max_c(BoundingBox2DCustom const& box) { return box[1]; }
  };

  using AdaptorCustom = OrthoTree::AdaptorGeneralBase<2, Point2DCustom, BoundingBox2DCustom, AdaptorBasicsCustom, float>;


  // Tailored Quadtree objects

  using QuadtreePointCustom = OrthoTree::OrthoTreePoint<2, Point2DCustom, BoundingBox2DCustom, AdaptorCustom, float>;
  using QuadtreeBoxCustom = OrthoTree::OrthoTreeBoundingBox<2, Point2DCustom, BoundingBox2DCustom, AdaptorCustom, float>;
```

## Benchmarks
Octree creation for 3 point sets using different placing strategy, and Cylindrical point set generation time:<br>
<div align="center" width="100%"><img src="https://github.com/attcs/Octree/blob/master/docs/octree_point_perf_3sets.png" align="center"></div>
<br><br>
Octree creation for 3 box sets using different placing strategy, and Cylindrical box set generation time:<br><br>
<div align="center" width="100%"><img src="https://github.com/attcs/Octree/blob/master/docs/octree_box_perf_3sets.png" align="center"></div>
<br><br>
Collision detection:<br><br>
<div align="center" width="100%"><img src="https://github.com/attcs/Octree/blob/master/docs/collisiondetection.png" align="center"></div>
<br>
<br>
*CPU: AMD Ryzen 5 5600X 6-Core @ 3.70GHz, CPU benchmark: 22146

