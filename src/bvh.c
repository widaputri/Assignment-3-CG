#include "bvh.h"
#include <string.h>
#include <stdio.h>
#include <assert.h>

// Comparison function for qsort
typedef struct {
    uint32_t axis;
    const Primitive* primitives;
} SortContext;

static SortContext sort_ctx;

static int compare_primitives(const void* a, const void* b) {
    uint32_t idx_a = *(const uint32_t*)a;
    uint32_t idx_b = *(const uint32_t*)b;

    Vec3 center_a = aabb_center(sort_ctx.primitives[idx_a].bounds);
    Vec3 center_b = aabb_center(sort_ctx.primitives[idx_b].bounds);

    float val_a = ((float*)&center_a)[sort_ctx.axis];
    float val_b = ((float*)&center_b)[sort_ctx.axis];

    if (val_a < val_b) return -1;
    if (val_a > val_b) return 1;
    return 0;
}

// Find best split using SAH
SplitCandidate bvh_find_best_split(const BVH* bvh, uint32_t* prim_indices,
                                   uint32_t start, uint32_t end) {
    SplitCandidate best = {FLT_MAX, 0, start + (end - start) / 2};
    const uint32_t num_bins = 12;

    for (uint32_t axis = 0; axis < 3; axis++) {
        // Compute bounds for this subset
        AABB bounds = aabb_empty();
        for (uint32_t i = start; i < end; i++) {
            bounds = aabb_union(bounds, bvh->primitives[prim_indices[i]].bounds);
        }

        float axis_min = ((float*)&bounds.min)[axis];
        float axis_max = ((float*)&bounds.max)[axis];

        if (axis_max - axis_min < 0.0001f) continue;

        // Binning
        typedef struct {
            AABB bounds;
            uint32_t count;
        } Bin;

        Bin bins[12] = {0};
        float bin_width = (axis_max - axis_min) / num_bins;

        // Fill bins
        for (uint32_t i = start; i < end; i++) {
            Vec3 center = aabb_center(bvh->primitives[prim_indices[i]].bounds);
            float pos = ((float*)&center)[axis];
            uint32_t bin_idx = (uint32_t)((pos - axis_min) / bin_width);
            if (bin_idx >= num_bins) bin_idx = num_bins - 1;

            bins[bin_idx].count++;
            if (bins[bin_idx].count == 1) {
                bins[bin_idx].bounds = bvh->primitives[prim_indices[i]].bounds;
            } else {
                bins[bin_idx].bounds = aabb_union(bins[bin_idx].bounds,
                                                  bvh->primitives[prim_indices[i]].bounds);
            }
        }

        // Sweep to find best split
        for (uint32_t split_bin = 1; split_bin < num_bins; split_bin++) {
            AABB left_bounds = aabb_empty();
            AABB right_bounds = aabb_empty();
            uint32_t left_count = 0;
            uint32_t right_count = 0;

            for (uint32_t i = 0; i < split_bin; i++) {
                if (bins[i].count > 0) {
                    left_bounds = aabb_union(left_bounds, bins[i].bounds);
                    left_count += bins[i].count;
                }
            }

            for (uint32_t i = split_bin; i < num_bins; i++) {
                if (bins[i].count > 0) {
                    right_bounds = aabb_union(right_bounds, bins[i].bounds);
                    right_count += bins[i].count;
                }
            }

            if (left_count == 0 || right_count == 0) continue;

            // Calculate SAH cost
            float traversal_cost = 1.0f;
            float intersect_cost = 1.0f;
            
            float left_area = aabb_surface_area(left_bounds);
            float right_area = aabb_surface_area(right_bounds);
            float parent_area = aabb_surface_area(bounds);
            
            float cost = traversal_cost + 
                        (left_count * left_area + right_count * right_area) / parent_area * intersect_cost;
            
            // Update best split if this is better
            if (cost < best.cost) {
                best.cost = cost;
                best.split_axis = axis;
                
                // Find the actual split position
                float split_pos = axis_min + split_bin * bin_width;
                
                // Partition primitives
                uint32_t left_idx = start;
                for (uint32_t i = start; i < end; i++) {
                    Vec3 center = aabb_center(bvh->primitives[prim_indices[i]].bounds);
                    float pos = ((float*)&center)[axis];
                    if (pos < split_pos) {
                        uint32_t temp = prim_indices[left_idx];
                        prim_indices[left_idx] = prim_indices[i];
                        prim_indices[i] = temp;
                        left_idx++;
                    }
                }
                
                best.split_pos = left_idx;
            }
        }
    }

    return best;
}

// Build BVH recursively
BVHNode* bvh_build_recursive(BVH* bvh, uint32_t* prim_indices,
                            uint32_t start, uint32_t end, uint32_t* node_idx) {
    // Allocate new node
    BVHNode* node = &bvh->nodes[(*node_idx)++];

    // Compute bounding box for all primitives in this node
    node->bounds = aabb_empty();
    for (uint32_t i = start; i < end; i++) {
        node->bounds = aabb_union(node->bounds, bvh->primitives[prim_indices[i]].bounds);
    }

    uint32_t prim_count = end - start;

    // Create leaf node if primitive count is small enough
    if (prim_count <= 2) {
        node->is_leaf = true;
        node->first_prim_idx = start;
        node->prim_count = prim_count;
        return node;
    }

    // Find best split using SAH
    SplitCandidate split = bvh_find_best_split(bvh, prim_indices, start, end);
    
    // Check if split is valid
    if (split.split_pos <= start || split.split_pos >= end) {
        // Fallback to median split using qsort if SAH failed
        uint32_t longest_axis = 0;
        Vec3 extent = vec3_sub(node->bounds.max, node->bounds.min);
        if (extent.y > extent.x && extent.y > extent.z) longest_axis = 1;
        else if (extent.z > extent.x) longest_axis = 2;
        
        sort_ctx.axis = longest_axis;
        sort_ctx.primitives = bvh->primitives;
        qsort(&prim_indices[start], prim_count, sizeof(uint32_t), compare_primitives);
        
        split.split_pos = start + prim_count / 2;
        
        // If still invalid, create leaf
        if (split.split_pos <= start || split.split_pos >= end) {
            node->is_leaf = true;
            node->first_prim_idx = start;
            node->prim_count = prim_count;
            return node;
        }
    }
    
    // Create internal node
    node->is_leaf = false;
    
    // Recursively build left and right subtrees
    node->left = bvh_build_recursive(bvh, prim_indices, start, split.split_pos, node_idx);
    node->right = bvh_build_recursive(bvh, prim_indices, split.split_pos, end, node_idx);

    return node;
}

// Create BVH
BVH* bvh_create(Primitive* primitives, uint32_t count) {
    BVH* bvh = (BVH*)calloc(1, sizeof(BVH));
    bvh->primitives = primitives;
    bvh->prim_count = count;

    // Allocate nodes (worst case: 2N-1 nodes)
    bvh->nodes = (BVHNode*)calloc(2 * count - 1, sizeof(BVHNode));
    bvh->indices = (uint32_t*)malloc(count * sizeof(uint32_t));

    // Initialize indices
    for (uint32_t i = 0; i < count; i++) {
        bvh->indices[i] = i;
    }

    // Build tree
    uint32_t node_idx = 0;
    bvh->root = bvh_build_recursive(bvh, bvh->indices, 0, count, &node_idx);
    bvh->node_count = node_idx;

    // Reorder primitives according to indices
    Primitive* reordered = (Primitive*)malloc(count * sizeof(Primitive));
    for (uint32_t i = 0; i < count; i++) {
        reordered[i] = primitives[bvh->indices[i]];
    }
    memcpy(primitives, reordered, count * sizeof(Primitive));
    free(reordered);

    return bvh;
}

// Destroy BVH
void bvh_destroy(BVH* bvh) {
    if (bvh) {
        free(bvh->nodes);
        free(bvh->indices);
        free(bvh);
    }
}

// BVH traversal (iterative for performance)
bool bvh_hit(const BVH* bvh, const Ray* ray, float t_min, float t_max,
             HitRecord* rec) {
    // Stack for iterative traversal
    BVHNode* stack[64];
    int stack_ptr = 0;

    bool hit_anything = false;
    float closest_so_far = t_max;

    // Start with root node
    if (bvh->root) {
        stack[stack_ptr++] = bvh->root;
    }

    // Traverse the BVH tree
    while (stack_ptr > 0) {
        // Pop node from stack
        BVHNode* node = stack[--stack_ptr];
        
        // Test AABB intersection
        if (!aabb_hit(&node->bounds, ray, t_min, closest_so_far)) {
            continue;
        }
        
        if (node->is_leaf) {
            // Test all primitives in this leaf node
            for (uint32_t i = 0; i < node->prim_count; i++) {
                uint32_t prim_idx = node->first_prim_idx + i;
                if (primitive_hit(&bvh->primitives[prim_idx], ray, t_min, closest_so_far, rec)) {
                    hit_anything = true;
                    closest_so_far = rec->t;
                }
            }
        } else {
            // Internal node - push children to stack
            // Push in order to visit closer child first (if possible)
            if (node->left) {
                stack[stack_ptr++] = node->left;
            }
            if (node->right) {
                stack[stack_ptr++] = node->right;
            }
        }
    }

    return hit_anything;
}