#include <vector>
#include <math.h>
#include <iostream>
#include <algorithm>
#include "pose_trees/Octree/octree.h"
#include "pose_trees/Octree/octree_container.h"

using namespace OrthoTree;

class TwinOcTrees
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
    TwinOcTrees()
    {
        c_tree_plus = RotationalTreePointND<6>();
        c_tree_minus = RotationalTreePointND<6>();
    }

    void create(vector<array<double, 6>> input_poses, int depth, int poses_in_leaf_node)
    {
        poses_plus = input_poses;
        poses_minus = input_poses;

        std::for_each(poses_minus.begin(), poses_minus.end(), [](array<double, 6> &n)
        {
            for (size_t i = 3; i < 6; i++)
            {
                if (n[i] < 180)
                {
                    n[i] = 180.0 + n[i];
                }
                else
                {
                    n[i] = n[i] - 180.0;
                }
            }
        });

        std::array<double, 6> inspection_space_min = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        std::array<double, 6> inspection_space_max = {100.0, 100.0, 100.0, 359.0, 359.0, 359.0};
        OrthoTree::BoundingBoxND<6> inspection_space;
        inspection_space.Min = inspection_space_min;
        inspection_space.Max = inspection_space_max;

        c_tree_plus.Create(c_tree_plus, poses_plus, depth, inspection_space, poses_in_leaf_node);

        c_tree_minus.Create(c_tree_minus, poses_minus, depth, inspection_space, poses_in_leaf_node);
    }

    std::pair<double,int> getNearestNeighbor(array<double, 6> pose)
    {
        return getNearestNeighbors(pose, 1)[0];
    }

    //Not very pretty and there is probaply a more performant version, but it works.
    //NNS of both Trees with concatenating of the two results
    vector<std::pair<double,int>> getNearestNeighbors(array<double, 6> search_point, int k)
    {
        vector<long unsigned int> neighbors_plus = c_tree_plus.GetNearestNeighbors(search_point, k, poses_plus);

        array<double, 6> search_point_revert = search_point;

        for (size_t i = 3; i < 6; i++)
        {

            if (search_point_revert[i] < 180)
            {
                search_point_revert[i] = 180.0 + search_point_revert[i];
            }
            else
            {
                search_point_revert[i] = search_point_revert[i] - 180.0;
            }
        }
        vector<long unsigned int> neighbors_minus = c_tree_minus.GetNearestNeighbors(search_point_revert, k, poses_minus);

        vector<std::pair<double,int>> results;

        vector<vector<long unsigned int>> neighbor_containers;
        neighbor_containers.push_back(neighbors_minus);
        for(int i = 0; i < k; i++)
        {
            double pose_distance;
            int neighbor;
            int neighbor_plus = neighbors_plus[i];
            double tranlation_distance_plus = euclidean_distance(search_point[0], search_point[1], search_point[2], poses_plus[neighbor_plus][0], poses_plus[neighbor_plus][1], poses_plus[neighbor_plus][2]);
            double rot_distance_plus = rotaional_distance(search_point[3], search_point[4], search_point[5], poses_plus[neighbor_plus][3], poses_plus[neighbor_plus][4], poses_plus[neighbor_plus][5]);
            double pose_distance_plus = pose_metric(tranlation_distance_plus, rot_distance_plus);
            std::pair<double,int>P = std::make_pair(pose_distance_plus, neighbor_plus);
            results.push_back(P);
        }
        for (size_t i = 0; i < neighbor_containers.size(); i++)
        {
            vector<long unsigned int> neighbor_container = neighbor_containers[i];
            for (size_t j = 0; j < k; j++)
            {
                int neighbor_minus = neighbor_container[j];
                double tranlation_distance_minus = euclidean_distance(search_point_revert[0], search_point_revert[1], search_point_revert[2], poses_minus[neighbor_minus][0], poses_minus[neighbor_minus][1], poses_minus[neighbor_minus][2]);
                double rot_distance_minus = rotaional_distance(search_point_revert[3], search_point_revert[4], search_point_revert[5], poses_minus[neighbor_minus][3], poses_minus[neighbor_minus][4], poses_minus[neighbor_minus][5]);
                double pose_distance_minus = pose_metric(tranlation_distance_minus, rot_distance_minus);
                bool allowed_to_add = true;
                for (size_t l = 0; l < results.size(); l++)
                {
                    auto pair = results[l];
                    if (pair.second == neighbor_minus)
                    {
                        allowed_to_add = false;
                        if (pair.first > pose_distance_minus )
                        {
                            pair.first = pose_distance_minus;
                        }
                    }
                }
                if (allowed_to_add)
                {
                    std::pair<double,int>P = std::make_pair(pose_distance_minus, neighbor_minus);
                    results.push_back(P);
                }
            }
        }
        std::sort(results.begin(),results.end());
        return results;
    }


private:
    RotationalTreePointND<6> c_tree_plus;
    RotationalTreePointND<6> c_tree_minus;
    vector<array<double, 6>> poses_plus;
    vector<array<double, 6>> poses_minus;
};
