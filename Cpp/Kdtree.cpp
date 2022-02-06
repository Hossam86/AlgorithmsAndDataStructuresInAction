#include <iostream>
#include <random>
#include <vector>
struct vec3
{
    float x;
    float y;
    float z;
};

bool operator==(vec3 &p1, vec3 &p2)
{
    return (p1.x == p2.x && p1.y == p2.y && p1.z == p2.z);
}

struct KdNode
{
    size_t level;
    vec3 point;
    KdNode *left;
    KdNode *right;
};

float get_point_key(vec3 p, size_t level, size_t k)
{
    return (&p.x)[level % k];
}

float get_node_key(KdNode *n, size_t k)
{
    return get_point_key(n->point, n->level, k);
}

// compare key of point and key of node
int compare(vec3 p, KdNode *n, size_t k)
{
    return get_point_key(p, n->level, k) - get_node_key(n, k);
}

// compue the distance between apoint and its projection on the split line passing through an node
float splitDistance(vec3 p, KdNode *n, size_t k)
{
    return std::abs(get_point_key(p, n->level, k) - get_node_key(n, k));
}

KdNode *
search(KdNode *root, KdNode *target)
{
    if (root == nullptr)
        return nullptr;
    else if (root->point == target->point)
        return root;
    else if (compare(root->point, target, 3) < 0)
        return search(root->left, target);
    else
        return search(root->right, target);
}

KdNode *
insert(KdNode *n, vec3 p, size_t level)
{
    if (n == nullptr)
    {
        KdNode *new_node = new KdNode{};
        new_node->level = level;
        new_node->left = nullptr;
        new_node->right = nullptr;
        return new_node;
    }

    else if (n->point == p)
    {
        return n;
    }

    else if (compare(p, n, 3) < 0)
    {
        n->left = insert(n->left, p, n->level + 1);
        return n;
    }
    else
    {
        n->right = insert(n->right, p, n->level + 1);
        return n;
    }
}

inline static void
swap(KdNode *x, KdNode *y)
{
    vec3 tmp = x->point;
    x->point = y->point;
    y->point = tmp;
}

KdNode *
_quickselect_median(KdNode *start, KdNode *end, int idx)
{
    if (end <= start)
        return nullptr;
    if (end == start + 1)
        return start;

    KdNode *md = start + (end - start) / 2;
    while (true)
    {
        float pivot = (&md->point.x)[idx];

        swap(md, end - 1);
        KdNode *store = start;
        for (auto p = start; p < end; p++)
        {
            if ((&p->point.x)[idx] < pivot)
            {
                if (p != store)
                    swap(p, store);
                store++;
            }
        }
        swap(store, end - 1);

        // median has duplicate values
        if ((&store->point.x)[idx] == (&md->point.x)[idx])
            return md;

        if (store > md)
            end = store;
        else
            start = store;
    }
}

static KdNode *
_construct_median_space_split(KdNode *node, uint32_t len, int i)
{
    if (len == 0)
        return nullptr;

    KdNode *n = _quickselect_median(node, node + len, i);
    if (n != nullptr)
    {
        i = (i + 1) % 3;
        n->left = _construct_median_space_split(node, n - node, i);
        n->right = _construct_median_space_split(n + 1, node + len - (n + 1), i);
    }
    return n;
}

struct Kdtree
{
    std::vector<KdNode> nodes;
    KdNode *root;
    size_t k;
};

void
Kdtree_free(Kdtree &self)
{
    self.nodes.clear();
}

Kdtree
kdtree_build(vec3 *vertices, size_t count, int axis)
{
    if (vertices == nullptr)
        return Kdtree{};

    Kdtree self{};
    for (size_t i = 0; i < count; ++i)
    {
        KdNode *tmp=new KdNode{};
        tmp->point = vertices[i];
        tmp->left = nullptr;
        tmp->right = nullptr;
        self.nodes.push_back(*tmp);
    }

    self.root = _construct_median_space_split(&(self.nodes[0]), count, axis);
    return self;
}

KdNode *
Kdtree_search(Kdtree tree, vec3 p)
{
    KdNode target;
    target.point = p;
    return search(tree.root, &target);
}
float squared_dist(vec3 p1, vec3 p2)
{
    return (p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y) + (p2.z - p1.z) * (p2.z - p1.z);
}
void nearest_neighbor(KdNode *root, const vec3 &target, KdNode *nearest_node, float *best_dist)
{
    if (root == nullptr)
        return;
    else
    {
        float dist = squared_dist(root->point, target);
        if (dist < *best_dist)
        {
            *best_dist = dist;
            *nearest_node = *root;
        }
        KdNode *close_branch;
        KdNode *far_branch;

        if (compare(target, root, 3) < 0)
        {
            close_branch = root->left;
            far_branch = root->right;
        }
        else
        {
            close_branch = root->right;
            far_branch = root->left;
        }

        nearest_neighbor(close_branch, target, nearest_node, best_dist);
        if (splitDistance(target, root, 3) < *best_dist)
            nearest_neighbor(far_branch, target, nearest_node, best_dist);
    }
}

int main(int argc, char const *argv[])
{

    std::vector<vec3> points;
    std::mt19937 rng;

    rng.seed(0);
    for (size_t level = 0; level < 3; ++level)
    {
        float r = (level + 1) * 10.0f;
        for (size_t i = 0; i < 5; i++)
        {
            std::uniform_int_distribution<> dis(0, 100);
            auto u = dis(rng) / 100.0f;
            auto v = dis(rng) / 100.0f;

            float theta = 2 * u * 3.14;
            float phi = 2 * v * 3.14;

            float rx = r * sin(phi) * cos(theta);
            float ry = r * sin(phi) * sin(theta);
            float rz = r * cos(phi);
            points.push_back(vec3{rx, ry, rz});
        }
    }
    auto query_point = vec3{0.0, 0.0, 0.0};
    auto tree = kdtree_build(&points[0], points.size(), 0);
    float best_dist = std::numeric_limits<float>::max();
    KdNode *nearest_node = new KdNode{};

    nearest_neighbor(tree.root, query_point, nearest_node, &best_dist);
    std::cout << squared_dist(query_point, nearest_node->point) << std::endl;
    return 0;
}
