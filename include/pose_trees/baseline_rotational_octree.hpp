#include <vector>
#include <math.h>
#include <iostream>
#include "pose_trees/Octree/octree.h"
#include "pose_trees/Octree/octree_container.h"

using namespace OrthoTree;

class BaselineRotationalOctree
{
private:
    double pose_metric(double tranlation_distance, double rot_distance)
    {
        return (tranlation_distance + rot_distance);
    }

    double euclidean_distance(double x1, double y1, double z1, double x2, double y2, double z2)
    {
        return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
    }

    double angle_distance(double query_angle, double angle2)
    {
        double distance = abs(query_angle - angle2);
        if (distance > 180)
        {
            if (query_angle < angle2)
            {
                return distance = abs((query_angle + 360) - angle2);
            }
            else
            {
                return distance = abs(query_angle - (angle2 + 360));
            }
        }
        return distance;
    }

    double rotaional_distance(double angleX1, double angleY1, double angleZ1, double angleX2, double angleY2, double angleZ2)
    {
        double diff_x = angle_distance(angleX1, angleX2);
        double diff_y = angle_distance(angleY1, angleY2);
        double diff_z = angle_distance(angleZ1, angleZ2);

        return ((diff_x + diff_y + diff_z) / 3);
    }

public:
    BaselineRotationalOctree()
    {
        array<array<double, 3>, 0> constexpr points;
        c_tree_ = OctreePoint();
    }

    void create(vector<array<double, 6>> poses, int depth, int poses_in_leaf_node)
    {
        depth_ = depth;
        nr_of_poses = poses.size();
        for (size_t i = 0; i < poses.size(); i++)
        {
            array<double,3> point = {poses[i][0], poses[i][1], poses[i][2]};
            points.insert(points.end(), point);
            array<double,3> ori = {poses[i][3], poses[i][4], poses[i][5]};
            orientations.insert(orientations.end(), ori);
        }


        std::array<double, 3> inspection_space_min = {0.0, 0.0, 0.0};
        std::array<double, 3> inspection_space_max = {100.0, 100.0, 100.0};
        OrthoTree::BoundingBoxND<3> inspection_space;
        inspection_space.Min = inspection_space_min;
        inspection_space.Max = inspection_space_max;
        c_tree_.Create(c_tree_, points, depth, inspection_space, poses_in_leaf_node);

        //Create a lookup vector for the K-overhead required.
        std::vector<std::pair<int, int>> table = {
            {0, 1},
            {100, 70},
            {1000, 310},
            {5000, 830},
            {10000, 1300},
            {20000, 1800},
            {50000, 2900},
            {100000, 5000}
        };

        for (size_t i = 0; i < table.size() - 1; i++)
        {
            int poseStart = table[i].first;
            int poseEnd = table[i +1].first;

            int kStart = table[i].second;
            int kEnd = table[i +1].second;

            for (int pose = poseStart; pose <= poseEnd; pose++)
            {
                double interpolationFactor = static_cast<double>(pose - poseStart) / (poseEnd - poseStart);
                int interpolatedValue = static_cast<int>(kStart + interpolationFactor * (kEnd - kStart));
                interpolatedArray[pose] = interpolatedValue;
            }
        }
    }

    vector<int> getNearestNeighbors(array<double, 6> pose, int k)
    {
        int k_with_overhead = interpolatedArray[nr_of_poses];

        array<double, 3> point = {pose[0], pose[1], pose[2]};
        vector<long unsigned int> nearestNeighbors = c_tree_.GetNearestNeighbors(point, k_with_overhead, points);
        vector<std::pair<double,int>> results;
        vector<int> results_ids;

        for(int i = 0; i < nearestNeighbors.size(); i++)
        {
            int neighbor = nearestNeighbors[i];
            double tranlation_distance = euclidean_distance(pose[0], pose[1], pose[2], points[neighbor][0], points[neighbor][1], points[neighbor][2]);
            double rot_distance = rotaional_distance(pose[3], pose[4], pose[5], orientations[neighbor][0], orientations[neighbor][1], orientations[neighbor][2]);
            double pose_distance = pose_metric(tranlation_distance, rot_distance);

            std::pair<double,int>P = std::make_pair(pose_distance, neighbor);
            results.push_back(P);
        }
        std::sort(results.begin(),results.end());

        for (size_t i = 0; i < k; i++)
        {
            results_ids.push_back(results[i].second);
        }

        return results_ids;
    }

private:
    OctreePoint c_tree_;
    int depth_;
    int interpolatedArray[100001];          //Biggest entry of the lookup vector +1
    int nr_of_poses;
    vector<array<double, 3>> orientations;
    vector<array<double, 3>> points;
};