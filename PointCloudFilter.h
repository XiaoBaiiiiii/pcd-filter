#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <cmath>
#include <algorithm>
#include <queue>
#include <string>
#include <sstream>
#include <pcl/io/ply_io.h>
#include <pcl/point_types.h>
#include <pcl/filters/extract_indices.h>

namespace pointcloudfilter {

template< class Real > 
struct OctreeNode;

template< class Real > 
struct Point {
    Point(size_t pt_index, Real pt_x, Real pt_y, Real pt_z) {
        index = pt_index;
        x = pt_x;
        y = pt_y;
        z = pt_z;}

    size_t index;
    Real x, y, z;
    double gaussian_curv;
    OctreeNode<Real>* node = nullptr;
    void Node(OctreeNode<Real>* n) {node = n;}
};

template< class Real > 
struct OctreeNode {

    Real centroid[3];
    Real back_left_bottom[3];
    Real front_right_top[3];
    Real half_size;
    OctreeNode* children[8];
    OctreeNode* father_node = nullptr;
    std::map<size_t, Point<Real>> points;
    bool leaf_node = true;

    OctreeNode(Real center[3], Real size ) {
        centroid[0] = center[0];
        centroid[1] = center[1];
        centroid[2] = center[2];

        half_size = size/2;

        back_left_bottom[0] = center[0] - size/2;
        back_left_bottom[1] = center[1] - size/2;
        back_left_bottom[2] = center[2] - size/2;
        front_right_top[0] = center[0] + size/2;
        front_right_top[1] = center[1] + size/2;
        front_right_top[2] = center[2] + size/2;
    }

    size_t PointNum() const {return points.size();}
    void fatherNode(OctreeNode* n) {father_node = n;}
    void createChildren(){
        // create children nodes
        Real center_1[3] = {centroid[0] + half_size/2, centroid[1] + half_size/2, centroid[2] + half_size/2};              
        children[0] = new OctreeNode(center_1, half_size);

        Real center_2[3] = {centroid[0] - half_size/2, centroid[1] + half_size/2, centroid[2] + half_size/2};             
        children[1] = new OctreeNode(center_2, half_size);

        Real center_3[3] = {centroid[0] + half_size/2, centroid[1] - half_size/2, centroid[2] + half_size/2};                 
        children[2] = new OctreeNode(center_3, half_size);
  
        Real center_4[3] = {centroid[0] - half_size/2, centroid[1] - half_size/2, centroid[2] + half_size/2};           
        children[3] = new OctreeNode(center_4, half_size);
     
        Real center_5[3] = {centroid[0] + half_size/2, centroid[1] + half_size/2, centroid[2] - half_size/2};           
        children[4] = new OctreeNode(center_5, half_size);

        Real center_6[3] = {centroid[0] - half_size/2, centroid[1] + half_size/2, centroid[2] - half_size/2};               
        children[5] = new OctreeNode(center_6, half_size);
      
        Real center_7[3] = {centroid[0] + half_size/2, centroid[1] - half_size/2, centroid[2] - half_size/2};               
        children[6] = new OctreeNode(center_7, half_size);
 
        Real center_8[3] = {centroid[0] - half_size/2, centroid[1] - half_size/2, centroid[2] - half_size/2};               
        children[7] = new OctreeNode(center_8, half_size);

        for (int i = 0; i < 8; i++) {
            children[i]->fatherNode(this);
        }

        leaf_node = false;
    }
    // if a given point coordinate in the bounding box
    bool InBoundingBox(const Point<Real>& point) {
        if (point.x < back_left_bottom[0] || 
            point.x > front_right_top[0] || 
            point.y < back_left_bottom[1] || 
            point.y > front_right_top[1] || 
            point.z < back_left_bottom[2] || 
            point.z > front_right_top[2]) {
            return false;
        } 
        return true;

    }
    void insertPoint(const Point<Real>& point) {points.insert(std::make_pair(point.index, point));}
    void deletePoint(Point<Real>& point) {
        auto& iter = points.find(point.index);
        if (iter != points.end()) {points.erase(iter);}
    }

};

template <class Real>
class Comp{
    public:
    bool operator()(const std::pair<double,Point<Real>>&n1,const std::pair<double,Point<Real>>&n2) const
    {
            if(n1.first == n2.first)
                    return n1.first > n2.first;
            return n1.first > n2.first;
    }
};


template< class Real >
class PointCloud {
public:
    virtual void getPoint(size_t index, Real* coords) = 0;
    virtual size_t size() = 0;
    // waiting for coding
    const std::vector<Point<Real>> allPoints(){
        std::vector<Point<Real>> all_points;
        size_t index = 0;
        Real* coords;
        while (index < size()) {
            getPoint(index,coords);
            Point<Real> pt(index,coords[0],coords[1],coords[2]);
            all_points.push_back(pt);
        }
        return all_points;
    }

    // virtual void deletePoint(size_t index) = 0;
    // PointCloud& deletePoints(std::vector<Points>& points) {
    //     PointCloud cloud_out;
    //     cloud_out = *this;
    //     for (Point& pt : points) {
    //         size_t i = pt.index;
    //         cloud_out.deletePoint(i);
    //     }
    //     return cloud_out;
    // }

};

template <class Real, typename PointCloudType>
class PointCloudFilter {
public:
    // PointCloudFilter(PointCloudType& pc)cloud_(pc){}
    ~PointCloudFilter(){}
    // input the gaussian curvature, resolution and cloud
    // output point indexs that shouldn't be used.

    std::set<size_t> planeDownSample(const double& gaussian_curvature, 
                        const double& resolution,
                        const PointCloudType& cloud_in){

        // Build octree
        OctreeNode<Real>* rootNode = buildOctree(cloud_in.allPoints());
        std::set<size_t> points_delete;
        for (const Point<Real>& pt: cloud_in.allPoints()) {
            std::vector<Point<Real>> neighbors;
            if (findKNearestNeighbors(8,pt,pt.node,neighbors)){
                double curv = calculateGaussianCurvature(pt, neighbors);
                if (curv < 0.01) {
                    for (auto& pt_neighbor : neighbors) {
                        double d = std::sqrt(std::pow(pt_neighbor.x - pt.x,2) + std::pow(pt_neighbor.y - pt.y,2) + std::pow(pt_neighbor.z - pt.z,2));
                        if (d < resolution) {
                            points_delete.insert(pt.index);
                            deletePoint(pt, pt.node);
                            break;
                        }
                    }
                } else {
                    continue;
                }
            } else {
                continue;
            }
        }
        return points_delete;
    }

    std::vector<size_t> planeDownSample(const double& gaussian_curvature, 
                        const double& resolution,
                        std::vector<Point<Real>>& points_in){

        // Build octree
        // point in the points_in store the leaf node in the octree
        OctreeNode<Real>* rootNode = buildOctree(points_in);

        std::vector<size_t> points_delete;

        for (auto& pt: points_in) {
            std::vector<Point<Real>> neighbors;
            if (findKNearestNeighbors(8,pt,pt.node,neighbors)){
                double curv = calculateGaussianCurvature(pt, neighbors);
                if (curv < 0.01) {
                    for (auto& pt_neighbor : neighbors) {
                        double d = std::sqrt(std::pow(pt_neighbor.x - pt.x,2) + std::pow(pt_neighbor.y - pt.y,2) + std::pow(pt_neighbor.z - pt.z,2));
                        if (d < resolution) {
                            points_delete.push_back(pt.index);
                            deletePoint(pt, pt.node);
                            break;
                        }
                    }
                } else {
                    continue;
                }
            } else {
                continue;
            }
        }
        return points_delete;
    }
    
    

private:

    OctreeNode<Real>* buildOctree(std::vector<Point<Real>>& points_in){
        // Calculate bounding box
        Real minX = points_in[0].x;
        Real maxX = points_in[0].x;
        Real minY = points_in[0].y;
        Real maxY = points_in[0].y;
        Real minZ = points_in[0].z;
        Real maxZ = points_in[0].z;

        for (auto& point : points_in)
        {
            minX = std::min(minX, point.x);
            maxX = std::max(maxX, point.x);
            minY = std::min(minY, point.y);
            maxY = std::max(maxY, point.y);
            minZ = std::min(minZ, point.z);
            maxZ = std::max(maxZ, point.z);
        }

        Real center[3] = {(minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2};

        Real size = std::max({maxX - minX, maxY - minY, maxZ - minZ});

        OctreeNode<Real>* rootNode = new OctreeNode<Real>(center, size);

        for (auto& point : points_in)
        {
            insertPoint(point, rootNode);
        }

        return rootNode;
    }

    // Build octree recursively
    void insertPoint(Point<Real>& point, OctreeNode<Real>* node, OctreeNode<Real>* father_node = nullptr) {
        // Check if point is in node's bounding box

        if (!node->InBoundingBox(point)){return;}

        if (node->PointNum() >= 8 && ((node->half_size) > 0.05)){
            node->insertPoint(point);

            if (node->leaf_node) {

                // create children nodes
                node->createChildren();
                
                for (auto& pt : node->points) {
                    for (int i = 0; i < 8; i++){
                        if (node->children[i]->InBoundingBox(pt.second)) {
                            insertPoint(pt.second, node->children[i], node);
                            break;
                        }
                    }
                }
            } else {
                for (int i = 0; i < 8; i++){
                    if (node->children[i]->InBoundingBox(point)) {
                        insertPoint(point, node->children[i],node);
                        break;
                    }
                }
            }


        } else {
            node->insertPoint(point);
            point.Node(node);


        }

    }

    // delete point from octree
    bool deletePoint(const Point<Real>& point,OctreeNode<Real>* node) {
        auto iter = node->points.find(point.index);
        if (iter != node->points.end()){
            node->points.erase(iter);
            if (node->father_node != nullptr){
                return deletePoint(point,node->father_node);
            }
        } else {
            std::cout<<"search point error"<<std::endl;
            return false;
        }
        return true;
    }

    /***
     * @param k number of the nearest neighbors
     * @param point point position used to calculate curvature
     * @param node point's node used to find neighbors
     * @param neighbors store the return value
     */
    // find k nearest neighbors
    bool findKNearestNeighbors(int k, const Point<Real>& point, OctreeNode<Real>* node, std::vector<Point<Real>>& neighbors) {
        if (node == nullptr) return false;
        if (node->PointNum()==0) return false;
        if (node->PointNum() > k) {

            std::priority_queue<std::pair<double, Point<Real>>,std::vector<std::pair<double, Point<Real>>>,Comp<Real>> distance_queue;

            for (auto iter = node->points.begin(); iter != node->points.end(); iter++) {
                if (iter->first == point.index) continue;
                double distance = static_cast<double>(std::sqrt(std::pow(point.x - iter->second.x, 2) + std::pow(point.y - iter->second.y, 2) + std::pow(point.z - iter->second.z, 2)));
                distance_queue.push(std::make_pair(distance,iter->second));
                if (distance_queue.size() == k+1) {
                    distance_queue.pop();
                }
            }
            if (distance_queue.size()!=k) return false;
            while (!distance_queue.empty()) {
                neighbors.push_back(distance_queue.top().second);
                distance_queue.pop();
            }
            return true;

        } else {
            return findKNearestNeighbors(k,point,node->father_node,neighbors);
        }

    }

    double calculateGaussianCurvature(const Point<Real>& pt, const std::vector<Point<Real>>& pt_vec) {
        const double eps = 1e-10; // 一个极小的值，用于判断分母是否为0
        double sum_uu = 0, sum_uv = 0, sum_vv = 0;
        double u, v, x, y, z;
        u = v = 0;

        for (const Point<Real>& pt_i : pt_vec) {
            x = pt_i.x - pt.x;
            y = pt_i.y - pt.y;
            z = pt_i.z - pt.z;

            u = std::atan2(y, x);
            v = std::acos(z / std::sqrt(x * x + y * y + z * z));

            double du = -y / (x * x + y * y);
            double dv = -(x * z) / (std::pow(x * x + y * y + z * z, 1.5));

            double duu = (2 * x * y) / std::pow(x * x + y * y, 2);
            double duv = (x * z) / (std::pow(x * x + y * y, 1.5) * std::sqrt(1 - z * z / (x * x + y * y)));
            double dvv = (x * x + y * y) * z / (std::pow(x * x + y * y + z * z, 2.5));

            sum_uu += duu * duu;
            sum_uv += duu * dv;
            sum_vv += dv * dv;
        }

        double k = (sum_uu * sum_vv - sum_uv * sum_uv) / (std::pow(sum_uu + sum_vv, 2) + eps);
        return k;
    }

    PointCloudType cloud_;
    OctreeNode<Real>* octree_;
    Real octree_min_resolution_ = 0.1;
};


std::vector<Point<float>> readPointCloud(const std::string& filename) {
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::io::loadPLYFile(filename, *cloud);
    std::vector<Point<float>> points;
    points.reserve(cloud->size());
    for (size_t i = 0; i < cloud->size(); i++) {
        const pcl::PointXYZ pt = cloud->points[i];
        points.emplace_back(i, pt.x, pt.y, pt.z);
    }
    return points;
}

void writePointCloud(const std::string& filename, pcl::PointCloud<pcl::PointXYZINormal>::Ptr cloud, const std::vector<size_t>& indices_to_remove) {

    pcl::PointIndicesPtr indices(new pcl::PointIndices);
    indices->indices.resize(indices_to_remove.size());

    for (size_t i = 0; i < indices_to_remove.size(); ++i) {
        indices->indices[i] = indices_to_remove[i];
    }

    pcl::PointCloud<pcl::PointXYZINormal>::Ptr filtered_cloud(new pcl::PointCloud<pcl::PointXYZINormal>);
    pcl::ExtractIndices<pcl::PointXYZINormal> extract;
    extract.setInputCloud(cloud);
    extract.setIndices(indices);
    extract.setNegative(true);
    extract.filter(*filtered_cloud);
    std::cout<<"number of points write out: "<< filtered_cloud->size() <<std::endl;

    pcl::io::savePLYFile(filename, *filtered_cloud);
}

bool pointCloudFilter(const std::string& file_in, 
                      const std::string& file_out, 
                      const float& resolution = 0.05, 
                      const int& k = 8, 
                      const double& gaussian_curv = 0.01) {

    std::vector<Point<float>> point_vec = readPointCloud(file_in);
    PointCloudFilter<float,pcl::PointCloud<pcl::PointXYZINormal>>* filter;
    std::vector<size_t> points_filtered_vec = filter->planeDownSample(gaussian_curv, resolution, point_vec);
    std::cout<<"number of points read in: "<<point_vec.size() <<std::endl;
    std::cout<<"number of points need to be filtered: "<< points_filtered_vec.size() <<std::endl;

    pcl::PointCloud<pcl::PointXYZINormal>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZINormal>);

    pcl::io::loadPLYFile(file_in, *cloud);
    writePointCloud(file_out, cloud, points_filtered_vec);
}

}// pointcloudfilter