#include <algorithm>
#include <cassert>
#include "BVH.hpp"

BVHAccel::BVHAccel(std::vector<Object *> p, int maxPrimsInNode,
                   SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod),
      primitives(std::move(p))
{
    if (primitives.empty())
        return;
    root = recursiveBuild(primitives);
}

BVHBuildNode *BVHAccel::recursiveBuild(std::vector<Object *> objects)
{
    BVHBuildNode *node = new BVHBuildNode();

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->first->getBounds());
    if (objects.size() == 1)
    {
        // Create leaf _BVHBuildNode_
        node->bounds = objects[0]->first->getBounds();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        return node;
    }
    else if (objects.size() == 2)
    {
        node->left = recursiveBuild(std::vector<Object *>{objects[0]});
        node->right = recursiveBuild(std::vector<Object *>{objects[1]});

        node->bounds = Union(node->left->bounds, node->right->bounds);
        return node;
    }
    else
    {
        Bounds3 centroidBounds;
        for (int i = 0; i < objects.size(); ++i)
            centroidBounds =
                Union(centroidBounds, objects[i]->first->getBounds().Centroid());
        int dim = centroidBounds.maxExtent();
        switch (dim)
        {
        case 0:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->first->getBounds().Centroid()[0] <
                       f2->first->getBounds().Centroid()[0];
            });
            break;
        case 1:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->first->getBounds().Centroid()[1] <
                       f2->first->getBounds().Centroid()[1];
            });
            break;
        case 2:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->first->getBounds().Centroid()[2] <
                       f2->first->getBounds().Centroid()[2];
            });
            break;
        }

        auto beginning = objects.begin();
        auto middling = objects.begin() + (objects.size() / 2);
        auto ending = objects.end();

        auto leftshapes = std::vector<Object *>(beginning, middling);
        auto rightshapes = std::vector<Object *>(middling, ending);

        assert(objects.size() == (leftshapes.size() + rightshapes.size()));

        node->left = recursiveBuild(leftshapes);
        node->right = recursiveBuild(rightshapes);

        node->bounds = Union(node->left->bounds, node->right->bounds);
    }

    return node;
}

bool BVHAccel::Intersect(const Ray &ray, double t0, double t1, HitRecord &hr) const
{
    if (!root)
        return false;
    bool isect = BVHAccel::getIntersection(root, ray, t0, t1, hr);
    return isect;
}

bool BVHAccel::getIntersection(BVHBuildNode *node, const Ray &ray, double t0, double t1, HitRecord &hr) const
{
    // TODO Traverse the BVH to find intersection
    bool hit = false;
    if (!node)
        return hit;
    if (node->isLeaf())
    {
        if (node->object->first->intersect(ray, t0, t1, hr))
        {
            hit = true;
            t1 = hr.t;
            hr.f = node->object->second;
            // if (hr.n.x() == 0.0 && hr.n.y() == 0.0 && hr.n.z() == 0.0)
            //     hr.n = SlVector3(0.0, 0.0, 1.0);
        }
    }
    else if (node->bounds.IntersectP(ray))
    {
        HitRecord hr1, hr2;
        bool isect1 = getIntersection(node->left, ray, t0, t1, hr1);
        if (isect1)
            t1 = hr1.t;
        bool isect2 = getIntersection(node->right, ray, t0, t1, hr2);
        if (!isect1 && !isect2)
        {
            hit = false;
        }
        else
        {
            hr = (isect1 && isect2) ? (hr1.t < hr2.t ? hr1 : hr2) : (isect1 ? hr1 : hr2);
            hit = true;
        }
    }
    return hit;
}