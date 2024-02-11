#include <chrono>
#include <random>
#include <Eigen/Geometry>

#include "../include/pose_trees/baseline_rotational_octree.hpp"
#include "../include/pose_trees/twin_octrees.hpp"
#include "../include/pose_trees/CSVWriter.h"

auto now = std::chrono::system_clock::now();
auto UTC = std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()).count();

std::string save_path = "";

double inner_product(array<double, 7> quat1, array<double, 7> quat2)
{
    double product = 0.0;
    for (size_t iDim = 3; iDim < 7; ++iDim)
    product += quat1[iDim] * quat2[iDim];

    return product;
}

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

array<double, 3> angle_to_rad(array<double, 3> pose)
{
    constexpr double pi = 3.14159265358979323846;
    array<double, 3> result;
    for (size_t i = 0; i < 3; i++)
    {
        result[i] = pose[i] * pi / 180;
    }
    return result;
}
Eigen::Quaternionf normalize_quaternion(Eigen::Quaternionf q)
{
    Eigen::Quaternionf normalized_q;
    auto norm = pow(q.x(),2) + pow(q.y(),2) + pow(q.z(),2) + pow(q.w(),2);
    auto s = pow(norm, (-0,5));
    normalized_q.x() = q.x() * s;
    normalized_q.y() = q.y() * s;
    normalized_q.z() = q.z() * s;
    normalized_q.w() = q.w() * s;

    return normalized_q;
}

vector<array<double, 7>> euler_to_quat(vector<array<double, 6>> poses)
{
    vector<array<double, 7>> quat_poses;
    for (size_t i = 0; i < poses.size(); i++)
    {
        array<double, 3> angle_in_rad = angle_to_rad({poses[i][3], poses[i][4], poses[i][5]});
        Eigen::Quaternionf q;
        q = Eigen::AngleAxisf(angle_in_rad[0], Eigen::Vector3f::UnitX()) * Eigen::AngleAxisf(angle_in_rad[1], Eigen::Vector3f::UnitY()) * Eigen::AngleAxisf(angle_in_rad[2], Eigen::Vector3f::UnitZ());
        Eigen::Quaternionf normalized_q = normalize_quaternion(q);
        quat_poses.push_back({poses[i][0], poses[i][1], poses[i][2], normalized_q.w(), normalized_q.x(), normalized_q.y(), normalized_q.z()});
    }
    return quat_poses;
}

void test_setup(int number_of_poses, int tree_depth, int number_of_neighbors_k, int poses_in_leaf_node, int number_of_test, CSVWriter& csv_file)
{

    double performance_listsearch;
    double performance_baseline;
    double performance_twin_tree_container;
    double performance_quaternion_tree;
    double performance_quat_listsearch;
    double performance_algorithm;
    double performance_algorithm_V2;
    double performance_6_dim_tree;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> trans_dist(1.0, 99.0);
    std::uniform_real_distribution<double> rot_dist(1.0, 358.0);

    std::uniform_real_distribution<double> search_rot_dist(2.0, 358.0);

    vector<array<double, 6>> poses;

    for (size_t i = 0; i < number_of_poses; i++)
    {
        array<double, 6> pose_3d = {trans_dist(mt), trans_dist(mt), trans_dist(mt), rot_dist(mt), rot_dist(mt), rot_dist(mt)};
        poses.push_back(pose_3d);
    }
    array<double, 6> search_point = {trans_dist(mt), trans_dist(mt), trans_dist(mt), search_rot_dist(mt), search_rot_dist(mt), search_rot_dist(mt)};

    // Quaternion setup
    vector<array<double, 7>> quat_poses = euler_to_quat(poses);
    array<double, 7> quat_search_point = euler_to_quat({search_point})[0];

    {
        //List search
        auto start = std::chrono::high_resolution_clock::now();                         // CLOCK START
        vector<std::pair<double,int>> results;
        for (size_t i = 0; i < number_of_poses; i++)
        {
            array<double,6> pose = poses[i];
            double tranlation_distance = euclidean_distance(search_point[0], search_point[1], search_point[2], pose[0], pose[1], pose[2]);
            double rot_distance = rotaional_distance(search_point[3], search_point[4], search_point[5], pose[3], pose[4], pose[5]);
            double pose_distance = pose_metric(tranlation_distance, rot_distance);

            std::pair<double,int>P = std::make_pair(pose_distance, i);
            results.push_back(P);
        }
        std::sort(results.begin(),results.end());
        auto stop = std::chrono::high_resolution_clock::now();                           // CLOCK END

        auto duration_micro = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        performance_listsearch = duration_micro.count();
    }

    {
        //Baseline_rotational_Octree
        auto tree = BaselineRotationalOctree();
        tree.create(poses, tree_depth, poses_in_leaf_node);

        auto start = std::chrono::high_resolution_clock::now();                             // CLOCK START
        vector<int> neighbors = tree.getNearestNeighbors(search_point, number_of_neighbors_k);
        auto stop = std::chrono::high_resolution_clock::now();                              // CLOCK END

        auto duration_micro = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        performance_baseline = duration_micro.count();
    }

    {
        //TwinTree Container
        auto tree = TwinOcTrees();
        tree.create(poses, tree_depth, poses_in_leaf_node);

        auto start = std::chrono::high_resolution_clock::now();                             // CLOCK START
        vector<std::pair<double,int>> neighbors = tree.getNearestNeighbors(search_point, number_of_neighbors_k);
        auto stop = std::chrono::high_resolution_clock::now();                              // CLOCK END

        auto duration_micro = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        performance_twin_tree_container = duration_micro.count();

    }

    {
        //Smallest Node Algorithm
        std::array<double, 6> inspection_space_min = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        std::array<double, 6> inspection_space_max = {100.0, 100.0, 100.0, 359.0, 359.0, 359.0};
        OrthoTree::BoundingBoxND<6> inspection_space;
        inspection_space.Min = inspection_space_min;
        inspection_space.Max = inspection_space_max;

        auto tree = PoseTreePointND<6>();
        tree.Create(tree, poses, tree_depth, inspection_space, poses_in_leaf_node);

        auto start = std::chrono::high_resolution_clock::now();                             // CLOCK START
        auto neighbors = tree.GetNearestNeighbors(search_point, number_of_neighbors_k, poses);
        auto stop = std::chrono::high_resolution_clock::now();                              // CLOCK END

        auto duration_micro = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        performance_algorithm = duration_micro.count();
    }

    {
        //Box Distance Algorithm
        std::array<double, 6> inspection_space_min = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        std::array<double, 6> inspection_space_max = {100.0, 100.0, 100.0, 359.0, 359.0, 359.0};
        OrthoTree::BoundingBoxND<6> inspection_space;
        inspection_space.Min = inspection_space_min;
        inspection_space.Max = inspection_space_max;

        auto tree = PoseTreePointV2ND<6>();
        tree.Create(tree, poses, tree_depth, inspection_space, poses_in_leaf_node);

        auto start = std::chrono::high_resolution_clock::now();                             // CLOCK START
        auto neighbors = tree.GetNearestNeighbors(search_point, number_of_neighbors_k, poses);
        auto stop = std::chrono::high_resolution_clock::now();                              // CLOCK END

        auto duration_micro = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        performance_algorithm_V2 = duration_micro.count();
    }

    {
        //Quat List search
        auto start = std::chrono::high_resolution_clock::now();                         // CLOCK START
        vector<std::pair<double,int>> results;
        for (size_t i = 0; i < number_of_poses; i++)
        {
            array<double,7> quat_pose = quat_poses[i];
            auto trans_distance = euclidean_distance(quat_search_point[0], quat_search_point[1], quat_search_point[2], quat_pose[0], quat_pose[1], quat_pose[2]);
            auto quat_product = inner_product(quat_search_point, quat_pose);
            auto quat_dis = 180 * (1 - pow(quat_product, 2));
            auto pose_distance = (trans_distance + quat_dis);

            std::pair<double,int>P = std::make_pair(pose_distance, i);
            results.push_back(P);
        }
        std::sort(results.begin(),results.end());
        auto stop = std::chrono::high_resolution_clock::now();                           // CLOCK END

        auto duration_micro = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        performance_quat_listsearch = duration_micro.count();
    }

    {
        //Quaternion Tree
        std::array<double, 7> inspection_space_min = {0.0, 0.0, 0.0, -1.0, -1.0, -1.0, -1.0};
        std::array<double, 7> inspection_space_max = {100.0, 100.0, 100.0, 1.0, 1.0, 1.0, 1.0};
        OrthoTree::BoundingBoxND<7> inspection_space;

        auto tree = QuaternionTreePointND<7>();
        tree.Create(tree, quat_poses, tree_depth, inspection_space, poses_in_leaf_node);

        auto start = std::chrono::high_resolution_clock::now();                             // CLOCK START
        auto neighbors = tree.GetNearestNeighbors(quat_search_point, number_of_neighbors_k, quat_poses);
        auto stop = std::chrono::high_resolution_clock::now();                              // CLOCK END

        auto duration_micro = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        performance_quaternion_tree = duration_micro.count();

    }

    {
        //Standard 6-Dimensional Tree

        std::array<double, 6> inspection_space_min = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        std::array<double, 6> inspection_space_max = {100.0, 100.0, 100.0, 359.0, 359.0, 359.0};
        OrthoTree::BoundingBoxND<6> inspection_space;
        inspection_space.Min = inspection_space_min;
        inspection_space.Max = inspection_space_max;

        auto tree = TreePointND<6>();
        tree.Create(tree, poses, tree_depth, inspection_space, poses_in_leaf_node);

        auto start = std::chrono::high_resolution_clock::now();                             // CLOCK START
        auto neighbors = tree.GetNearestNeighbors(search_point, number_of_neighbors_k, poses);
        auto stop = std::chrono::high_resolution_clock::now();                              // CLOCK END

        auto duration_micro = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        performance_6_dim_tree = duration_micro.count();
    }

    csv_file.newRow() << performance_listsearch << performance_baseline << performance_twin_tree_container << performance_algorithm << performance_algorithm_V2 << performance_quat_listsearch << performance_quaternion_tree << performance_6_dim_tree << number_of_poses << tree_depth << number_of_neighbors_k << poses_in_leaf_node << number_of_test;
}

void do_histograms()
{
    //Parameters as vectors to permutate many combinations for mutliple Histograms
    vector<int> number_of_poses = {30000};
    vector<int> tree_depth = {9};
    vector<int> number_of_neighbors_k = {3};
    vector<int> poses_in_leaf_node = {21};
    int number_of_tests = 100;

    std::cout << "-----Test Setup-----" << std::endl;

    double blah = 0.0;
    for (size_t i = 0; i < number_of_poses.size(); i++)
    {
        for (size_t j = 0; j < tree_depth.size(); j++)
        {
            for (size_t k = 0; k < number_of_neighbors_k.size(); k++)
            {
                for (size_t l = 0; l < poses_in_leaf_node.size(); l++)
                {
                    CSVWriter csv;
                    csv.newRow() << "list_search" << "baseline" << "twin_tree" << "algorithm" << "algorithmV2" << "quaternion_listsearch" << "quaternion" << "6_dimensional_tree" << "number_of_poses" << "tree_depth" << "number_of_neighbors_k" << "poses_in_leaf_node" << "number_of_tests";
                    for (size_t m = 0; m < number_of_tests; m++)
                    {
                        double max_tests = number_of_tests * poses_in_leaf_node.size() * number_of_neighbors_k.size() * tree_depth.size() * number_of_poses.size();
                        std::cout<< "\x1B[2J\x1B[H"; // Clear screen
                        std::cout << "Test: " << blah << " of " << max_tests << " that's " << static_cast<double>(blah) / max_tests * 100 << "%" << std::endl;
                        std::cout << "Test with parameters: " << std::endl;
                        std::cout << "Number of poses: " << number_of_poses[i] << std::endl;
                        std::cout << "Tree depth: " << tree_depth[j] << std::endl;
                        std::cout << "Number of neighbors k: " << number_of_neighbors_k[k] << std::endl;
                        std::cout << "Poses in leaf node: " << poses_in_leaf_node[l] << std::endl;
                        std::cout << "Number of test: " << m << " of " <<  number_of_tests << std::endl;
                        test_setup(number_of_poses[i], tree_depth[j], number_of_neighbors_k[k], poses_in_leaf_node[l], m, csv);
                        blah++;
                    }
                    csv.writeToFile(save_path + "histogram_" + std::to_string(number_of_poses[i]) + "_" + std::to_string(tree_depth[j]) + "_" + std::to_string(number_of_neighbors_k[k]) + "_" + std::to_string(poses_in_leaf_node[l]) + "_at_time_"+ std::to_string(UTC) + ".csv");
                }
            }
        }
    }
}

void poses_overtime()
{
    int start_number_of_poses = 1000;                //500
    int end_number_of_poses = 100000;               //100.000
    int step_size = 1000;
    int tree_depth = 5;
    int number_of_neighbors_k = 3;
    int poses_in_leaf_node = 21;
    int number_of_tests = 10;

    std::cout << "-----Test Setup-----" << std::endl;
    CSVWriter csv;
    csv.newRow() << "list_search" << "baseline" << "twin_tree" << "algorithm" << "algorithmV2" << "quaternion_listsearch" << "quaternion" << "6_dimensional_tree" << "number_of_poses" << "tree_depth" << "number_of_neighbors_k" << "poses_in_leaf_node" << "number_of_tests";
    for (size_t i = start_number_of_poses; i <= end_number_of_poses; i+=step_size)
    {

        for (size_t j = 0; j < number_of_tests; j++)
        {
            std::cout<< "\x1B[2J\x1B[H"; // Clear screen
            std::cout << "Number of poses: " << i << " out of " << end_number_of_poses << std::endl;
            std::cout << "Test: " << j << " of " << number_of_tests << std::endl;
            test_setup(i, tree_depth, number_of_neighbors_k, poses_in_leaf_node, j, csv);
        }
        csv.newRow() << "---" << "---" << "---" << "---" << "---" << "---";
    }
    csv.writeToFile(save_path + "poses_overtime_" + std::to_string(start_number_of_poses) + "_" + std::to_string(end_number_of_poses) + "_" + std::to_string(step_size) + "_at_timeDELETE_"+ std::to_string(UTC) + ".csv");
}

void max_leaf_overtime()
{
    int number_of_poses = 30000;
    int tree_depth = 9;
    int number_of_neighbors_k = 3;
    int start_poses_in_leaf_node = 2;           //2
    int end_poses_in_leaf_node = 50;            //50
    int step_size = 1;
    int number_of_tests = 10;

    std::cout << "-----Test Setup-----" << std::endl;
    CSVWriter csv;
    csv.newRow() << "list_search" << "baseline" << "twin_tree" << "algorithm" << "algorithmV2" << "quaternion_listsearch" << "quaternion" << "6_dimensional_tree" << "number_of_poses" << "tree_depth" << "number_of_neighbors_k" << "poses_in_leaf_node" << "number_of_tests";
    for (size_t i = start_poses_in_leaf_node; i <= end_poses_in_leaf_node; i+=step_size)
    {

        for (size_t j = 0; j < number_of_tests; j++)
        {
            std::cout<< "\x1B[2J\x1B[H"; // Clear screen
            std::cout << "Number of poses in leaf node: " << i << " out of " << end_poses_in_leaf_node << std::endl;
            std::cout << "Test: " << j << " of " << number_of_tests << std::endl;
            test_setup(number_of_poses, tree_depth, number_of_neighbors_k, i, j, csv);
        }
        csv.newRow() << "---" << "---" << "---" << "---" << "---" << "---";
    }
    csv.writeToFile(save_path + "max_leaf_overtime_" + std::to_string(start_poses_in_leaf_node) + "_" + std::to_string(end_poses_in_leaf_node) + "_" + std::to_string(step_size) + "_at_time_"+ std::to_string(UTC) + ".csv");
}

void k_overtime()
{
    int number_of_poses = 30000;
    int tree_depth = 9;
    int start_number_of_neighbors_k = 1;            //1
    int end_number_of_neighbors_k = 70;             //70
    int step_size = 1;
    int poses_in_leaf_node = 21;
    int number_of_tests = 10;

    std::cout << "-----Test Setup-----" << std::endl;
    CSVWriter csv;
    csv.newRow() << "list_search" << "baseline" << "twin_tree" << "algorithm" << "algorithmV2" << "quaternion_listsearch" << "quaternion" << "6_dimensional_tree" << "number_of_poses" << "tree_depth" << "number_of_neighbors_k" << "poses_in_leaf_node" << "number_of_tests";
    for (size_t i = start_number_of_neighbors_k; i <= end_number_of_neighbors_k; i+=step_size)
    {

        for (size_t j = 0; j < number_of_tests; j++)
        {
            std::cout<< "\x1B[2J\x1B[H"; // Clear screen
            std::cout << "Number of neighbors k: " << i << " out of " << end_number_of_neighbors_k << std::endl;
            std::cout << "Test: " << j << " of " << number_of_tests << std::endl;
            test_setup(number_of_poses, tree_depth, i, poses_in_leaf_node, j, csv);
        }
        csv.newRow() << "---" << "---" << "---" << "---" << "---" << "---";
    }
    csv.writeToFile(save_path + "k_overtime_" + std::to_string(start_number_of_neighbors_k) + "_" + std::to_string(end_number_of_neighbors_k) + "_" + std::to_string(step_size) + "_at_time_"+ std::to_string(UTC) + ".csv");
}

void tree_depth_overtime()
{
    int number_of_poses = 30000;
    int start_tree_depth = 2;               //2
    int end_tree_depth = 9;                 //9
    int number_of_neighbors_k = 3;
    int step_size = 1;
    int poses_in_leaf_node = 21;
    int number_of_tests = 10;

    std::cout << "-----Test Setup-----" << std::endl;
    CSVWriter csv;
    csv.newRow() << "list_search" << "baseline" << "twin_tree" << "algorithm" << "algorithmV2" << "quaternion_listsearch" << "quaternion" << "6_dimensional_tree" << "number_of_poses" << "tree_depth" << "number_of_neighbors_k" << "poses_in_leaf_node" << "number_of_tests";
    for (size_t i = start_tree_depth; i <= end_tree_depth; i+=step_size)
    {

        for (size_t j = 0; j < number_of_tests; j++)
        {
            std::cout<< "\x1B[2J\x1B[H"; // Clear screen
            std::cout << "Tree depth: " << i << " out of " << end_tree_depth << std::endl;
            std::cout << "Test: " << j << " of " << number_of_tests << std::endl;
            test_setup(number_of_poses, i, number_of_neighbors_k, poses_in_leaf_node, j, csv);
        }
        csv.newRow() << "---" << "---" << "---" << "---" << "---" << "---";
    }
    csv.writeToFile(save_path + "depth_overtime_" + std::to_string(start_tree_depth) + "_" + std::to_string(end_tree_depth) + "_" + std::to_string(step_size) + "_at_time_"+ std::to_string(UTC) + ".csv");
}

int main()
{
    do_histograms();
    poses_overtime();
    max_leaf_overtime();
    k_overtime();
    tree_depth_overtime();
}