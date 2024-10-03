#include "wrapper.h"

#include <cstdlib>
#include <time.h>

#include <stdio.h>
#include <memory.h>
#include <cstdlib>

typedef struct KDNode KDNode;
typedef struct Bounds Bounds;

struct Bounds {
    float mins[3];
    float maxs[3];

    Bounds() {
        this->mins[0] = -(__FLT_MAX__);
        this->mins[1] = -(__FLT_MAX__);
        this->mins[2] = -(__FLT_MAX__);

        this->maxs[0] = __FLT_MAX__;
        this->maxs[1] = __FLT_MAX__;
        this->maxs[2] = __FLT_MAX__;
    }

    float dist2(float* point) {
        float totalDist2 = 0;

        for(int i = 0; i < 3; i++) {
            float diff = 0;

            if(point[i] < this->mins[i]) {
                diff = point[i] - this->mins[i];
            } else if(point[i] > this->maxs[i]) {
                diff = point[i] - this->maxs[i];
            }

            totalDist2 = totalDist2 + diff * diff;
        }

        return totalDist2;
    }

    void cut(Bounds* result, float* point, int comp, bool less) {
        result->mins[0] = this->mins[0];
        result->mins[1] = this->mins[1];
        result->mins[2] = this->mins[2];

        result->maxs[0] = this->maxs[0];
        result->maxs[1] = this->maxs[1];
        result->maxs[2] = this->maxs[2];

        
        if(less) {
            result->maxs[comp] = point[comp];
        } else {
            result->mins[comp] = point[comp];
        }
    }
};

struct KDNode {
    float point[3];
    KDNode* left;
    KDNode* right;
    int comp;
    Bounds bounds;

    void init(float* point, int comp) {
        this->point[0] = point[0];
        this->point[1] = point[1];
        this->point[2] = point[2];

        this->comp = comp;

        this->left = NULL;
        this->right = NULL;
    }

    float dist2(float* point) {
        float totalDist2 = 0;

        for(int i = 0; i < 3; i++) {
            float diff = this->point[i] - point[i];
            totalDist2 = totalDist2 + diff * diff;
        }

        return totalDist2;
    }

    void destroy() {
        if(this->left != NULL) {
            this->left->destroy();
            delete this->left;
        }

        if(this->right != NULL) {
            this->right->destroy();
            delete this->right;
        }
    }
};

class KDTree {
public:
    float m_minDist2;
    float* m_points;
    size_t m_size;
    size_t m_pointNum;
    size_t m_pointCount;
    KDNode* m_nodes;
    size_t m_nodeCount;

private:
    KDNode* allocNode(float* point, int comp) {
        KDNode* result = &m_nodes[m_nodeCount];
        result->init(point, comp);
        m_nodeCount++;
        return result;
    }

    void insertNode(KDNode* current, float* point, int comp) {
        if(m_nodeCount == 0) {
            allocNode(point, 0);
            return; 
        }

        int nextComp = (comp + 1) % 3;
        if(current->point[comp] > point[comp]) {
            if(current->left == NULL) {
                current->left = allocNode(point, nextComp);
                current->bounds.cut(&current->left->bounds, current->point, comp, true);
            } else {
                insertNode(current->left, point, nextComp);
            }
        } else {
            if(current->right == NULL) {
                current->right = allocNode(point, nextComp);
                current->bounds.cut(&current->right->bounds, current->point, comp, false);
            } else {
                insertNode(current->right, point, nextComp);
            }
        }
    }

    bool isTooClose(KDNode* current, float* point) {
        if(m_nodeCount == 0) return false;

        float dist2 = current->dist2(point);

        if(dist2 < m_minDist2) {
            return true;
        }

        if(current->left != NULL && current->left->bounds.dist2(point) < m_minDist2) {
            if(isTooClose(current->left, point)) {
                return true;
            }
        }

        if(current->right != NULL && current->right->bounds.dist2(point) < m_minDist2) {
            if(isTooClose(current->right, point)) {
                return true;
            }
        }

        return false;
    }

public:
    KDTree(float* points, size_t pointNum, float minDist2) {
        m_points = points;
        m_minDist2 = minDist2;
        m_pointNum = pointNum;
        m_pointCount = 0;

        m_nodes = (KDNode*)malloc(sizeof(KDNode) * (pointNum + 1));
        m_nodeCount = 0;
    }

    bool addPoint(float* point) {
        if(isTooClose(m_nodes, point)) {
            return false;
        }

        insertNode(m_nodes, point, 0);

        m_points[3*m_pointCount + 0] = point[0];
        m_points[3*m_pointCount + 1] = point[1];
        m_points[3*m_pointCount + 2] = point[2];

        m_pointCount++;

        return true;
    }
    
    void destroy() {
        free(m_nodes);
    }
};

void generate_seperated_points_extern(float* out, int element_num, float bound0, float bound1, float bound2, float min_atom_separation, int print_progress) {
    srand(time(NULL));

    KDTree* tree = new KDTree(out, element_num, min_atom_separation * min_atom_separation);

    float tempPos[3];
    
    for(int i = 0; i < element_num; i++) {
        tempPos[0] = bound0 * (((float)rand()) / ((float)RAND_MAX));
        tempPos[1] = bound1 * (((float)rand()) / ((float)RAND_MAX));
        tempPos[2] = bound2 * (((float)rand()) / ((float)RAND_MAX));

        if(!tree->addPoint(tempPos))
            i--;
        
        if(print_progress == 1 && i % 10000 == 0) printf("\rGenerated %d/%d positions", i, element_num);
    }

    if(print_progress == 1)
        printf("\n");
}