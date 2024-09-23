#include <chrono>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <cmath>
#include <cassert>

#include <glm/glm.hpp>
#include <raylib.h>

using namespace glm;
using std::vector;
using std::unordered_map;
using std::unordered_set;
// using std::unique_ptr;

#define __ITERATION_STRATEGY_SPATIAL_COHERENCE

std::vector<glm::vec2> collisions;
std::mutex entityMutex;
std::mutex partitioningMutex;
std::mutex already_processed_ents_mutex;

// #define pl(x) std::cout << ":" << __LINE__ << " " << __FUNCTION__ << "() " #x " " << x << "\n";
#define pl(x) 
#define apl(x) std::cout << ":" << __LINE__ << " " << __FUNCTION__ << "() " #x " " << x << "\n";

template<typename Type> bool noneEqual(const std::pair<Type&, Type&>& pair1, const std::pair<Type&, Type&>& pair2) {
    return !(pair1.first  == pair2.first || pair1.first  == pair2.second || 
             pair1.second == pair2.first || pair1.second == pair2.second);
}
bool aprox_equal(float x, float y, float eps = 0.0001f){
    return glm::abs(x-y) < abs(eps*x);
}
#define assert_apeq(x,y) do {\
    if(! aprox_equal(x,y)) {\
        pl(x)\
        pl(y)\
    }\
    assert(aprox_equal(x,y));\
} while(0)

#define TEST(fun) do {\
    std::cout << "entering" #fun << "\n";\
    fun;\
    std::cout << "leaving" #fun << "\n";\
}while(0)

// Do you really need more than that?
enum class ShapeType { 
    None, 
    Point, 
    Circle, 
    Square 
};

struct Entity {
    int id;
    vec2 position = vec2(0);
    vec2 velocity = vec2(0);
    float angularVelocity = 0.0f; // Added for rotation
    float rotation = 0.0f; // Rotation angle
    ShapeType shape = ShapeType::None; //TODO union
    union {
        struct{ // Circle
            float radius;
        };
        struct{ // Square
            float width = 0;
            float height = 0;
        };
    };
    float mass = 1;
    Entity* next = nullptr;
    int mip_stored = 0; // velocity AND size(rad/h+w) change mip, and to find entity after they cahnged we need this
    int idx_stored = 0; // Reserved
    std::function<void(Entity&, Entity&)> onCollision;

    float getCharacteristicSize(float deltaTime){
        // pl(radius)
        // pl((velocity * deltaTime).x)
        // pl((velocity * deltaTime).y)
        // pl(length(velocity * deltaTime))
        return (radius + length(velocity * deltaTime))*2.f; //TODO
    }
};

struct CollisionEvent {
    Entity* entity1 = nullptr;
    Entity* entity2 = nullptr;
    float timeUntilCollision = NAN;
};

struct SpatialCell {
    Entity* head = nullptr;
    std::mutex mutex;

    SpatialCell() = default;
    ~SpatialCell() = default;

    // Disable for mutex
    SpatialCell(const SpatialCell& other) {
        this->head = other.head;
        // this->mutex = other.mutex; TODO?
    }
    SpatialCell& operator=(const SpatialCell&) = delete;
    SpatialCell(SpatialCell&&) = delete;
    SpatialCell& operator=(SpatialCell&&) = delete;
};

static long int collisionCounter = 0;

class GridPartitioning {
public:
    ivec2 gridSize; //cells_x; cells_y
    float cellSize;
    vector<SpatialCell> cells;
    
    GridPartitioning(int gridWidth, int gridHeight, float _cellSize) 
        : gridSize(gridWidth, gridHeight), cellSize(_cellSize), cells(gridWidth * gridHeight) {}
    ~GridPartitioning() = default;

    int getCellIndex(const vec2& position) const {
        ivec2 xy = clamp(ivec2(position / vec2(cellSize)), ivec2(0), gridSize - 1);
        return xy.x + xy.y * gridSize.x;
    }

    void insertEntity(Entity* entity) {
        pl(entity->position.x)
        pl(entity->position.y)
        int index = getCellIndex(entity->position);
        pl(index)
        entity->idx_stored = index;
        entity->next = cells[index].head;
        cells[index].head = entity;
    }
    void insertEntityMT(Entity* entity) {
        int index = getCellIndex(entity->position);
        SpatialCell& cell = cells[index];
        std::lock_guard<std::mutex> lock(cell.mutex);  // Lock per cell for MT
        entity->idx_stored = index;
        entity->next = cell.head;
        cell.head = entity;
    }

    void removeEntity(Entity* entity) {
        int index = entity->idx_stored;
        Entity** current = &cells[index].head;
        while (*current != nullptr) {
            if (*current == entity) {
                *current = entity->next;
                entity->next = nullptr;
                return;
            }
            current = &(*current)->next;
        }
    }
    void removeEntityMT(Entity* entity) {
        int index = entity->idx_stored;
        SpatialCell& cell = cells[index];
        std::lock_guard<std::mutex> lock(cell.mutex);  // Lock per cell
        Entity** current = &cell.head;
        while (*current != nullptr) {
            if (*current == entity) {
                *current = entity->next;
                entity->next = nullptr;
                return;
            }
            current = &(*current)->next;
        }
    }

    // Tries to "cache" cell
    void moveEntity(Entity* entity) {
        int old_idx = entity->idx_stored;
        int new_idx = getCellIndex(entity->position);
        if(old_idx != new_idx){
            removeEntity(entity);
            insertEntity(entity);
        }
    }
    void moveEntityMT(Entity* entity) {
        int old_idx = entity->idx_stored;
        int new_idx = getCellIndex(entity->position);
        if (old_idx != new_idx) {
            {
                // std::lock_guard<std::mutex> lock(cells[old_idx].mutex);  // Lock old cell
                removeEntityMT(entity);
            }
            {
                // std::lock_guard<std::mutex> lock(cells[new_idx].mutex);  // Lock new cell
                insertEntityMT(entity);
            }
        }
    }

};

// None of the partitioning functions change entity physical properties
// This works effectively like a tree but little more predictable and limited
class MipPartitioning {
public:
    ivec2 baseGridSize;
    int baseCellSize;
    vector<GridPartitioning> grids;
    int mipLevels;
    
    // Base is one with highest dencity
    MipPartitioning(int baseGridWidth, int baseGridHeight, float _baseCellSize)
        : baseGridSize(baseGridWidth, baseGridHeight), baseCellSize(_baseCellSize) {
        mipLevels = glm::log2(float(max(baseGridHeight, baseGridWidth))) + 1;
        int currentWidth = baseGridWidth;
        int currentHeight = baseGridHeight;
        float currentCellSize = _baseCellSize;

        for (int mip = 0; mip < mipLevels; mip++) {
            grids.emplace_back(currentWidth, currentHeight, currentCellSize);
            currentWidth /= 2;
            currentHeight /= 2;
            currentCellSize *= 2;
        }
    }
    ~MipPartitioning() = default;

    int getMipLevel(float characteristicSize) {
        float ratio = characteristicSize / baseCellSize;
        int mipmap = clamp(int(glm::log2(ratio)),0, mipLevels-1);
        return mipmap;
        // return mipLevels-1;
    }

    // newCharSize should probably be ~ size + sqrt(dist moved on desired tick)
    void insertEntity(Entity* entity, float newCharSize) {
        int mip = getMipLevel(newCharSize); //TODO
        if(mip != 0) {pl(entity->radius)}
        if(mip != 0) {pl(length(entity->velocity))}
        // assert(mip == 0);
        
        pl(mip)
        assert(entity);
        entity->mip_stored = mip;
        grids[mip].insertEntity(entity);
    }
    void insertEntityMT(Entity* entity, float newCharSize) {
        int mip = getMipLevel(newCharSize); //TODO
        if(mip != 0) {pl(entity->radius)}
        if(mip != 0) {pl(length(entity->velocity))}
        // assert(mip == 0);
        
        pl(mip)
        assert(entity);
        entity->mip_stored = mip;
        grids[mip].insertEntityMT(entity);
    }

    void removeEntity(Entity* entity) {
        int mip = entity->mip_stored;
        grids[mip].removeEntity(entity);
    }
    void removeEntityMT(Entity* entity) {
        int mip = entity->mip_stored;
        grids[mip].removeEntityMT(entity);
    }

    void moveEntityAcrossMips(Entity* entity, float newCharSize){
        // Surely need to remove / insert. Different mips cannot share same cells
        removeEntity(entity);
        insertEntity(entity, newCharSize);
    }
    void moveEntityAcrossMipsMT(Entity* entity, float newCharSize){
        // Surely need to remove / insert. Different mips cannot share same cells
        removeEntityMT(entity);
        insertEntityMT(entity, newCharSize);
    }
    // In the same mip level
    void moveEntityAcrossCells(Entity* entity){
        int mip = entity->mip_stored;
        grids[mip].moveEntity(entity);
    }
    void moveEntityAcrossCellsMT(Entity* entity){
        int mip = entity->mip_stored;
        grids[mip].moveEntityMT(entity);
    }

    // Broad phase recalculation. Like remove+insert, but "cached"
    void moveEntity(Entity* entity, float newCharSize) {
        int new_mip = getMipLevel(newCharSize); //TODO
        if(entity->mip_stored == new_mip){
            moveEntityAcrossCells(entity);
        } else {
            moveEntityAcrossMips(entity, newCharSize);
            entity->mip_stored = new_mip;
        }
    }
    void moveEntityMT(Entity* entity, float newCharSize) {
        int new_mip = getMipLevel(newCharSize); //TODO
        if(entity->mip_stored == new_mip){
            moveEntityAcrossCellsMT(entity);
        } else {
            moveEntityAcrossMipsMT(entity, newCharSize);
            entity->mip_stored = new_mip;
        }
    }

    void destroy() {}
};

struct pair_hash {
    template <typename T1, typename T2>
    std::size_t operator()(const std::pair<T1, T2>& pair) const {
        auto hash1 = std::hash<T1>{}(pair.first);
        auto hash2 = std::hash<T2>{}(pair.second);
        // Combine two hashes symmetrically
        return hash1 + hash2;
    }
};

struct pair_equal {
    bool operator()(const std::pair<Entity*, Entity*>& lhs, const std::pair<Entity*, Entity*>& rhs) const {
        return ((lhs.first == rhs.first)  && (lhs.second == rhs.second)) || 
               ((lhs.first == rhs.second) && (lhs.second == rhs.first));
    }
};

struct ECS {
    std::unordered_map<int, Entity*> entities;
    //NEVER interract with entities in partitioning directly. This will lead to either crash, INF loop or memory leak
    MipPartitioning partitioning;

    ECS(int gridWidth, int gridHeight, float cellSize)
        : partitioning(gridWidth, gridHeight, cellSize) {}
    ~ECS() = default;
    
    void destroy() {}

    void insertEntity(int id, vec2 position, vec2 velocity, ShapeType shape, float size, float mass) {
        Entity* entity = new Entity{id, position, velocity, 0,0, shape, {0}, mass, nullptr};
        if (shape == ShapeType::Circle) {
            entity->radius = size;
        } else if (shape == ShapeType::Square) {
            entity->width = size;
            entity->height = size;
        }
        entities[id] = entity;
        partitioning.insertEntity(entity, entity->getCharacteristicSize(0));
    }

    // void updateEntityPosition(int id, float dtime) {
    //     Entity* entity = entities[id];
    //     partitioning.moveEntity(entity, entity->getCharacteristicSize(dtime));
    // }

    void processSimulationPrecise(float deltaTime) {
        float usedTime = 0.0f;
        std::unordered_set<std::pair<Entity*, Entity*>, pair_hash, pair_equal> already_processed_ents = {};
        
        // Continue finding and processing collisions until we've exhausted the delta time
        // this is like general physics engine approach but precise and painfully slow
        while (usedTime < deltaTime) {
            // pl("START ITER")
            CollisionEvent nextCollision = {0};
            bool collides_someday = findNextCollisionWithAll(deltaTime - usedTime, &nextCollision, already_processed_ents);
            bool collides_this_tick = (collides_someday) && (nextCollision.timeUntilCollision < (deltaTime - usedTime));
                pl(deltaTime)
                pl(usedTime)
                pl(collides_someday)
                pl(collides_this_tick)
                pl(nextCollision.timeUntilCollision)

            if(collides_this_tick){
                // Move all entities by the time until the next collision
                float collisionTime = min(nextCollision.timeUntilCollision, deltaTime-usedTime);

                if(collisionTime > 0){
                    moveAllEntitiesByTime(collisionTime);
                    usedTime += collisionTime;
                }

                // Handle the next collision
                if (nextCollision.entity1 && nextCollision.entity2) { //TODO redundant
                    vec2 delta = nextCollision.entity2->position - nextCollision.entity1->position;
                    resolveCollision(*nextCollision.entity1, *nextCollision.entity2, delta);
                    vec2 v = (
                        nextCollision.entity1->position * nextCollision.entity2->radius + 
                        nextCollision.entity2->position * nextCollision.entity1->radius
                        )/ (nextCollision.entity1->radius + nextCollision.entity2->radius);
                    collisions.push_back({v.x, v.y});
                }
            } else {
                moveAllEntitiesByTime(deltaTime - usedTime);
                break;
            }

            // Accumulate the used time
            // pl("END ITER")
        }
    }

    void processSimulationFastSingleThread(float dTime) {
        std::unordered_set<std::pair<Entity*, Entity*>, pair_hash, pair_equal> already_processed_ents = {};
        // for(int mip=0; mip<partitioning.mipLevels; mip++){
            #ifdef __ITERATION_STRATEGY_SPATIAL_COHERENCE
            for(int mip=0; mip<partitioning.mipLevels; mip++){
                for (int xx = 0; xx < partitioning.grids[mip].gridSize.x; xx++) {
                for (int yy = 0; yy < partitioning.grids[mip].gridSize.y; yy++) {
                    int idx = xx + yy * partitioning.grids[mip].gridSize.x;
                    Entity* entity = partitioning.grids[mip].cells[idx].head;
                    // Continue finding and processing collisions until we've exhausted the delta time
                    // this is like general physics engine approach but precise and painfully slow
                    while(entity) {
            #else
                    for (auto [id, entity]: entities) {
            #endif
                        // Continue finding and processing collisions until we've exhausted the delta time
                        // this is like general physics engine approach but precise and painfully slow
                        float usedTime = 0.0f;
                        while (usedTime < dTime) {
                            // pl("START ITER")
                            CollisionEvent nextCollision = {0};
                            bool collides_someday = findNextCollisionWithNearby(entity, &nextCollision, already_processed_ents);
                            bool collides_this_tick = (collides_someday) && (nextCollision.timeUntilCollision < (dTime - usedTime));
                                pl(deltaTime)
                                pl(usedTime)
                                pl(collides_someday)
                                pl(collides_this_tick)
                                pl(nextCollision.timeUntilCollision)

                            if(collides_this_tick){
                                // Move all entities by the time until the next collision
                                float collisionTime = min(nextCollision.timeUntilCollision, dTime-usedTime);

                                // if(collisionTime > 0){
                                    // moveAllEntitiesByTime(collisionTime);
                                    nextCollision.entity1->position += nextCollision.entity1->velocity * dTime;
                                    partitioning.moveEntity(nextCollision.entity1, nextCollision.entity1->getCharacteristicSize(dTime));
                                    
                                    nextCollision.entity2->position += nextCollision.entity2->velocity * dTime;
                                    partitioning.moveEntity(nextCollision.entity2, nextCollision.entity2->getCharacteristicSize(dTime));

                                    usedTime += collisionTime;

                                    // Handle the next collision
                                    if (nextCollision.entity1 && nextCollision.entity2) { //TODO redundant
                                        vec2 delta = nextCollision.entity2->position - nextCollision.entity1->position;
                                        resolveCollision(*nextCollision.entity1, *nextCollision.entity2, delta);
                                        vec2 v = (nextCollision.entity1->position + nextCollision.entity2->position)/2.f;
                                        collisions.push_back({v.x, v.y});
                                    }
                                // }
                            } else {
                                entity->position += entity->velocity * (dTime - usedTime);
                                partitioning.moveEntity(entity, entity->getCharacteristicSize(dTime - usedTime));
                                break;
                            }
                        }
            #ifdef __ITERATION_STRATEGY_SPATIAL_COHERENCE
                    entity = entity->next;
                    }
                }}
            }
            #else
            }
            #endif
    }
    
    
    void processSimulationFastThread(int startX, int endX, int startY, int endY, int mip, float dTime) {
        std::unordered_set<std::pair<Entity*, Entity*>, pair_hash, pair_equal> already_processed_ents;
        
        for (int xx = startX; xx < endX && xx < partitioning.grids[mip].gridSize.x; xx++) {
        for (int yy = startY; yy < endY && yy < partitioning.grids[mip].gridSize.y; yy++) {
            int idx = xx + yy * partitioning.grids[mip].gridSize.x;
            SpatialCell& cell = partitioning.grids[mip].cells[idx];
            {
                // std::lock_guard<std::mutex> lock(cell.mutex);  // Lock per cell
                Entity* entity = cell.head;
                while (entity) {
                    float usedTime = 0.0f;
                    while (usedTime < dTime) {
                        CollisionEvent nextCollision = {0};

                        bool collides_someday = findNextCollisionWithNearby(entity, &nextCollision, already_processed_ents);
                        bool collides_this_tick = (collides_someday) && (nextCollision.timeUntilCollision < (dTime - usedTime));

                        if (collides_this_tick) {
                            float collisionTime = min(nextCollision.timeUntilCollision, dTime - usedTime);

                            // Move entities safely (lock them during updates)
                            entityMutex.lock();
                            nextCollision.entity1->position += nextCollision.entity1->velocity * dTime;
                            partitioning.moveEntityMT(nextCollision.entity1, nextCollision.entity1->getCharacteristicSize(dTime));

                            nextCollision.entity2->position += nextCollision.entity2->velocity * dTime;
                            partitioning.moveEntityMT(nextCollision.entity2, nextCollision.entity2->getCharacteristicSize(dTime));
                            entityMutex.unlock();

                            usedTime += collisionTime;
                        } else {
                            entity->position += entity->velocity * (dTime - usedTime);
                            partitioning.moveEntityMT(entity, entity->getCharacteristicSize(dTime - usedTime));
                            break;
                        }
                    }
                    entity = entity->next;  // Proceed to the next entity in the cell
                }
            }
        }}
    }

    void moveAllEntitiesByTime(float dTime) {
        for (auto& [id, entity] : entities) {
            entity->position += entity->velocity * dTime;
            partitioning.moveEntity(entity, entity->getCharacteristicSize(dTime));
        }
    }

    // Resolve collision between two entities
    void resolveCollision(Entity& entity1, Entity& entity2, const vec2& delta) {
        vec2 normal = normalize(delta);
        vec2 relativeVelocity = entity2.velocity - entity1.velocity;
        float separatingVelocity = dot(relativeVelocity, normal);
        // if (!std::isnan(separatingVelocity)) std::cout << "separatingVelocity " << separatingVelocity << "\n";

        if (separatingVelocity < 0) {
            float impulseStrength = (2 * separatingVelocity) / (entity1.mass + entity2.mass);
            vec2 impulse = impulseStrength * normal;
            entity1.velocity += impulse * entity2.mass;
            entity2.velocity -= impulse * entity1.mass;
        }

        if (entity1.onCollision) entity1.onCollision(entity1, entity2);
        if (entity2.onCollision) entity2.onCollision(entity2, entity1);
    }

    bool findNextCollisionWithAll(float remainingTime, CollisionEvent* event, std::unordered_set<std::pair<Entity*, Entity*>, pair_hash, pair_equal>& already_processed_ents) {
        *event = {0};
        float closest_collision_time = +INFINITY;
        bool any_collision_found = false;
        
        // Check every pair of nearby entities and calculate when they'll collide
        // to find nearby ones, grab all the mip(map)s in 3x3 area. This IS slow.
        for(int mip=0; mip<partitioning.mipLevels; mip++){
            for (int xx = 0; xx < partitioning.grids[mip].gridSize.x; xx++) {
            for (int yy = 0; yy < partitioning.grids[mip].gridSize.y; yy++) {
                ivec2 xxyy = ivec2(xx, yy);

                int idx = xxyy.x + xxyy.y * partitioning.grids[mip].gridSize.x;
                Entity* entity = partitioning.grids[mip].cells[idx].head;

                while (entity) {
                    // Check every nearby cell in every bigger-cell mipmap
                    for(int _mip=0; _mip<partitioning.mipLevels; _mip++){
                        // Recalculated because different mip levels will have different id's for same point
                        // Calculating it directly from xxyy leads to presicion and update delay problems
                        // ivec2 low  = clamp((xxyy)/(1 << _mip)-1, ivec2(0), partitioning.grids[_mip].gridSize - 1);
                        // ivec2 high = clamp((xxyy)/(1 << _mip)+1, ivec2(0), partitioning.grids[_mip].gridSize - 1);

                        vec2 center = entity->position;
                        float cell_size = partitioning.grids[_mip].cellSize; 
                        ivec2 low  = clamp(ivec2((center-cell_size) / cell_size), ivec2(0), partitioning.grids[_mip].gridSize - 1);
                        ivec2 high = clamp(ivec2((center+cell_size) / cell_size), ivec2(0), partitioning.grids[_mip].gridSize - 1);
                        
                        // low = ivec2(0);
                        // high = partitioning.grids[_mip].gridSize-1;

                        for (int _x = low.x; _x <= high.x; _x++) {
                        for (int _y = low.y; _y <= high.y; _y++) {
                            int idx = _x + _y * partitioning.grids[_mip].gridSize.x;
                            assert(_mip < partitioning.grids.size());
                            Entity* other = partitioning.grids[_mip].cells[idx].head;
                            pl("")
                            pl(other)
                            if(other) pl(other->id)

                            while (other) {
                                // Calculate time until collision between entity and other
                                //TODO move checks around
                                if(other != entity){ //reduntant check
                                    // pl(entity)
                                    // pl(other)
                                    float timeUntilCollision;
                                    pl(entity->id)
                                    pl(other->id)
                                    bool collides = calculateCollisionTime(*entity, *other, &timeUntilCollision);
                                    // pl(collides)

                                    if(collides){
                                        pl("FOUND COLLISION")
                                        pl(timeUntilCollision)
                                        if ((timeUntilCollision < closest_collision_time)) {
                                            if(! already_processed_ents.contains({entity, other})){
                                                any_collision_found = true;
                                        // pl(timeUntilCollision);
                                                closest_collision_time = timeUntilCollision;
                                                *event = {entity, other, timeUntilCollision};
                                                pl(event->entity1)
                                                pl(event->entity2)
                                                already_processed_ents.insert({entity, other});
                                            }
                                        }
                                    }
                                }
                                other = other->next;
                            }
                        }}
                    }
                    entity = entity->next;
                }
            }
            }
        }
        collisionCounter += any_collision_found;
        return any_collision_found;
    }

    bool findNextCollisionWithNearby(Entity* target_entity, CollisionEvent* event, std::unordered_set<std::pair<Entity*, Entity*>, pair_hash, pair_equal>& already_processed_ents) {
        *event = {0};
        float closest_collision_time = +INFINITY;
        bool any_collision_found = false;
        
        // Check every pair of nearby entities and calculate when they'll collide
        // to find nearby ones, grab all the mip(map)s in 3x3 area. This IS slow.

        // Check every nearby cell in every bigger-cell mipmap
        for(int mip=0; mip<partitioning.mipLevels; mip++){
            vec2 center = target_entity->position;
            float cell_size = partitioning.grids[mip].cellSize; 
            ivec2 low  = clamp(ivec2((center-cell_size) / cell_size), ivec2(0), partitioning.grids[mip].gridSize - 1);
            ivec2 high = clamp(ivec2((center+cell_size) / cell_size), ivec2(0), partitioning.grids[mip].gridSize - 1);
            
            // low = ivec2(0);
            // high = partitioning.grids[_mip].gridSize-1;

            for (int _x = low.x; _x <= high.x; _x++) {
            for (int _y = low.y; _y <= high.y; _y++) {
                int idx = _x + _y * partitioning.grids[mip].gridSize.x;
                assert(mip < partitioning.grids.size());
                Entity* other = partitioning.grids[mip].cells[idx].head;
                pl("")
                pl(other)
                if(other) pl(other->id)

                while (other) {
                    // Calculate time until collision between entity and other
                    //TODO move checks around
                    if(other != target_entity){ //reduntant check
                        // pl(entity)
                        // pl(other)
                        float timeUntilCollision;
                        pl(entity->id)
                        pl(other->id)
                        bool collides = calculateCollisionTime(*target_entity, *other, &timeUntilCollision);
                        // pl(collides)

                        if(collides){
                            pl("FOUND COLLISION")
                            pl(timeUntilCollision)
                            if ((timeUntilCollision < closest_collision_time)) {
                                {
                                    std::lock_guard<std::mutex> lock(already_processed_ents_mutex);
                                
                                    if(! already_processed_ents.contains({target_entity, other})){
                                        any_collision_found = true;
                                        // pl(timeUntilCollision);
                                        closest_collision_time = timeUntilCollision;
                                        *event = {target_entity, other, timeUntilCollision};
                                        pl(event->entity1)
                                        pl(event->entity2)
                                        already_processed_ents.insert({target_entity, other});
                                    }
                                }
                            }
                        }
                    }
                    other = other->next;
                }
            }}
        }
        collisionCounter += any_collision_found;
        return any_collision_found;
    }
    
    bool circleCircleIntersectionTime(vec2 pos1, vec2 vel1, float radius1, 
                                       vec2 pos2, vec2 vel2, float radius2,
                                       float* timeUntilCollision) {
        // Relative position and velocity
        vec2 p_rel = pos2 - pos1;
        vec2 v_rel = vel2 - vel1;

        if(length(p_rel) <= radius1 + radius2) {
            *timeUntilCollision = 0;
            return true;
        }

        // pl(p_rel.x)
        // pl(v_rel.x)

        float r_sum = radius1 + radius2;  // Combined radii
        float r_sum_sq = r_sum * r_sum;   // Square of combined radii

        // Coefficients for the quadratic equation
        float a = dot(v_rel, v_rel);  // ||v_rel||^2
        float b = 2.0f * dot(p_rel, v_rel);  // 2 * dot(p_rel, v_rel)
        float c = dot(p_rel, p_rel) - r_sum_sq;  // ||p_rel||^2 - (r1 + r2)^2

        // Solve the quadratic equation: a * t^2 + b * t + c = 0
        float discriminant = b * b - 4 * a * c;

        if (discriminant < 0 || a == 0) {
            // No real intersection or no relative velocity
            return false;  // Return negative time to indicate no intersection
        }

        // Two possible times (t1 and t2)
        float sqrtDiscriminant = sqrt(discriminant);
        float t1 = (-b - sqrtDiscriminant) / (2.0f * a);
        float t2 = (-b + sqrtDiscriminant) / (2.0f * a);

        // We're interested in the first positive time
        if (t1 >= 0) {
            *timeUntilCollision = t1;
            // pl(t1)
            return true;
        }
        if (t2 >= 0) {
            *timeUntilCollision = t2;
            // pl(t2)
            return true;
        }

        // If both times are negative, there's no future collision
        return false;
    }

    bool squareSquareIntersectionTime(Entity& square1, Entity& square2, float* timeUntilCollision) {
        // AABB collision detection between rotated squares
        vec2 delta = square2.position - square1.position;
        if (abs(delta.x) < (square1.width / 2 + square2.width / 2) &&
            abs(delta.y) < (square1.height / 2 + square2.height / 2)) {
            *timeUntilCollision = 0; // Collision already occurred
            return true;
        }
        return false;
    }

    // Handle point collision with shapes (point-circle, point-square)
    bool pointShapeCollision(Entity& point, Entity& shape, float* timeUntilCollision) {
        if (shape.shape == ShapeType::Circle) {
            return circleCircleIntersectionTime(point.position, point.velocity, 0,
                                                shape.position, shape.velocity, shape.radius,
                                                timeUntilCollision);
        } else if (shape.shape == ShapeType::Square) {
            // Treat the point as a tiny square with zero size
            vec2 delta = shape.position - point.position;
            if (abs(delta.x) < shape.width / 2 && abs(delta.y) < shape.height / 2) {
                *timeUntilCollision = 0;
                return true;
            }
        }
        return false;
    }

    // Calculate time until collision between two entities
    bool calculateCollisionTime(Entity& entity1, Entity& entity2, float* timeUntilCollision) {
        vec2 relativeVelocity = entity2.velocity - entity1.velocity;
        vec2 relativePosition = entity2.position - entity1.position;

        // pl((int)entity1.shape)
        // pl((int)entity2.shape)
        // Calculate time to collision based on the type of shapes
        if (entity1.shape == ShapeType::Circle && entity2.shape == ShapeType::Circle) {
            return circleCircleIntersectionTime(entity1.position, entity1.velocity, entity1.radius,
                                                entity2.position, entity2.velocity, entity2.radius,
                                                timeUntilCollision);
        } else if (entity1.shape == ShapeType::Square && entity2.shape == ShapeType::Square) {
        // Square-Square collision logic
        return squareSquareIntersectionTime(entity1, entity2, timeUntilCollision);
        } else if (entity1.shape == ShapeType::Point || entity2.shape == ShapeType::Point) {
            // Point-Square or Point-Circle collision logic
            return pointShapeCollision(entity1, entity2, timeUntilCollision);
        }
        //TODO other
        
        // If no valid collision, return a negative value
        assert(false);
        return false;
    }

    // Handle the collision between two entities (similar to your previous implementation)
    // void handleCollision(Entity& entity1, Entity& entity2) {
    //     // Collision handling logic from previous code...
    //     vec2 delta = entity2.position - entity1.position;
    //     float distance = length(delta);

    //     if (entity1.shape == ShapeType::Circle && entity2.shape == ShapeType::Circle) {
    //         if (distance < entity1.radius + entity2.radius) {
    //             resolveCollision(entity1, entity2, delta, distance);
    //         }
    //     } else if ((entity1.shape == ShapeType::Point && entity2.shape == ShapeType::Circle) ||
    //                (entity1.shape == ShapeType::Circle && entity2.shape == ShapeType::Point)) {
    //         float circleRadius = entity1.shape == ShapeType::Circle ? entity1.radius : entity2.radius;
    //         if (distance < circleRadius) {
    //             resolveCollision(entity1, entity2, delta, distance);
    //         }
    //     } else if (entity1.shape == ShapeType::Square && entity2.shape == ShapeType::Square) {
    //         if (abs(delta.x) < (entity1.width / 2 + entity2.width / 2) &&
    //             abs(delta.y) < (entity1.height / 2 + entity2.height / 2)) {
    //             resolveCollision(entity1, entity2, delta, distance);
    //         }
    //     }
    // }

    // Too big deltaTimes dramatically slow everything down
    void tick_simulation_precise(float deltaTime) {
        // This is neccessary due to how entities are processed. Instead of steps, everything is analytically solved. 
        // When entities are too fast, processing kernel will not notice it on current mip level, so we move it on mip level with bigger cells
        for (auto& [id, entity] : entities) {
            // float max_step = length(entity->velocity * deltaTime);

            // int old_mip = partitioning.getMipLevel();
            // int old_mip = partitioning.getMipLevel(max_step);
            partitioning.moveEntity(entity, entity->getCharacteristicSize(deltaTime));
            // if()
                           
        }
        processSimulationPrecise(deltaTime);
    }

    void tick_simulation_fast(float deltaTime) {
        // This is neccessary due to how entities are processed. Instead of steps, everything is analytically solved. 
        // When entities are too fast, processing kernel will not notice it on current mip level, so we move it on mip level with bigger cells
        for (auto& [id, entity] : entities) {
            // float max_step = length(entity->velocity * deltaTime);

            // int old_mip = partitioning.getMipLevel();
            // int old_mip = partitioning.getMipLevel(max_step);
            partitioning.moveEntity(entity, entity->getCharacteristicSize(deltaTime));
            // if()
                           
        }
        processSimulationFastSingleThread(deltaTime);
    }

    void printEntityPositions() {
        for (const auto& [id, entity] : entities) {
            std::cout << "Entity " << id << ": Position (" << entity->position.x << ", " << entity->position.y << ")\n";
        }
    }
    void printEntityVelocities() {
        for (const auto& [id, entity] : entities) {
            std::cout << "Entity " << id << ": Velocity (" << entity->velocity.x << ", " << entity->velocity.y << ")\n";
        }
    }
    // full?
    void printEntityFulls() {
        for (const auto& [id, entity] : entities) {
            std::cout << "Entity " << id 
                << ": Position (" << entity->position.x << ", " << entity->position.y << ")\n"
                << ": Velocity (" << entity->velocity.x << ", " << entity->velocity.y << ")\n"
                << ": mip " << entity->mip_stored << "\n"
                << ": idx " << entity->idx_stored << "\n";
        }
    }

    void printEntityProblems() {
        int empty = 0;
        for(int mip=0; mip<partitioning.mipLevels; mip++){
            for (int xx = 0; xx < partitioning.grids[mip].gridSize.x; xx++) {
            for (int yy = 0; yy < partitioning.grids[mip].gridSize.y; yy++) {
                ivec2 xxyy = ivec2(xx, yy);
                int idx = xxyy.x + xxyy.y * partitioning.grids[mip].gridSize.x;
                Entity* entity = partitioning.grids[mip].cells[idx].head;
                int ctr = 0;
                while (entity) {
                    ctr++;
                    entity = entity->next;
                }
                if (ctr != 0){
                    std::cout << "ctr " << ctr << ": mip " << mip << "\n"
                    << ": idx " << idx << "\n";;
                }
                else empty++;
            }}
            std::cout << "empty " << empty << ": mip " << mip << "\n";
        }
    }
};

void benchmark(ECS& ecs, int numIterations) {
    auto startTime = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < numIterations; i++) {
        ecs.tick_simulation_precise(0.016f);
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = endTime - startTime;

    std::cout << "Benchmark: " << numIterations << " iterations took "
              << duration.count() << " seconds.\n";
    std::cout << (duration.count() / numIterations) * 1000 << " ms / iter\n";
}   

void test_no_collision_movement() {
    ECS ecs(10, 10, 50.0f);

    ecs.insertEntity(1, vec2(0.0f, 0.0f), vec2(1.0f, 0.0f), ShapeType::Circle, 5.0f, 1.0f);
    ecs.tick_simulation_precise(1.0f);  // Simulate 1 second of movement
    // ecs.printEntityPositions();

    Entity* entity = ecs.entities[1];
    assert_apeq(entity->position.x, +1.0f);
}

void test_basic_collision() {
    ECS ecs(1, 1, +INFINITY);

    ecs.insertEntity(1, vec2( 0.0000f, 0.0f), vec2(+1.0f, 0.0f), ShapeType::Circle, 1.0f, .001f);
    ecs.insertEntity(2, vec2(10.0001f, 0.0f), vec2(-1.0f, 0.0f), ShapeType::Circle, 1.0f, 10000000.0f);

    ecs.tick_simulation_precise(100000000.0);
    
    // ecs.printEntityPositions();

    Entity* entity1 = ecs.entities[1];
    Entity* entity2 = ecs.entities[2];
    // std::cout << "entity1->velocity.x " << entity1->velocity.x << "\n";
    // std::cout << "entity2->velocity.x " << entity2->velocity.x << "\n";
    assert_apeq(entity1->velocity.x, -3);
    assert_apeq(entity2->velocity.x, -1);

    assert_apeq(entity1->position.x, -3e+08);
    assert_apeq(entity2->position.x, -1e+08);
}

void test_mass_impact_on_velocity() {
    ECS ecs(10, 10, 50.0f);

    ecs.insertEntity(1, vec2(0.0f, 0.0f), vec2(1.0f, 0.0f), ShapeType::Circle, 5.0f, 1.0f);
    ecs.insertEntity(2, vec2(10.0f, 0.0f), vec2(-1.0f, 0.0f), ShapeType::Circle, 5.0f, 10.0f);

    ecs.tick_simulation_precise(5.0f); // Simulate 5 seconds, expect a collision halfway
    // ecs.runCollisions();
    // ecs.printEntityPositions();

    Entity* entity1 = ecs.entities[1];
    Entity* entity2 = ecs.entities[2];

    // It is not obvious how this should look like so just general check
    if (entity1->velocity.x < 0.0f && std::abs(entity2->velocity.x) < std::abs(entity1->velocity.x)) {
        // std::cout << "Test passed: Mass correctly affected velocity.\n";
    } else {
        std::cout << "Test failed: Mass did not impact velocity as expected.\n";
    }
}

void test_circle_circle_collision_time() {
    ECS ecs(10, 10, 5.0f);

    ecs.insertEntity(1, vec2(0.0f, 0.0f), vec2(2.0f, 0.0f), ShapeType::Circle, 5.0f, 1.0f);
    ecs.insertEntity(2, vec2(20.0f, 0.0f), vec2(-2.0f, 0.0f), ShapeType::Circle, 5.0f, 1.0f);

    ecs.tick_simulation_precise(2.5f);
    // Expect the entities to collide after 2.5 seconds
    Entity* entity1 = ecs.entities[1];
    Entity* entity2 = ecs.entities[2];
    assert_apeq(entity1->position.x, 5.0f);
    assert_apeq(entity2->position.x, 15.0f);

    // ecs.printEntityPositions();
}

void test_entity_removal() {
    ECS ecs(10, 10, 50.0f);

    ecs.insertEntity(1, vec2(0.0f, 0.0f), vec2(1.0f, 0.0f), ShapeType::Circle, 5.0f, 1.0f);
    ecs.insertEntity(2, vec2(10.0f, 0.0f), vec2(-1.0f, 0.0f), ShapeType::Circle, 5.0f, 1.0f);

    ecs.tick_simulation_precise(0.0f);  // Simulate no movement

    ecs.partitioning.removeEntity(ecs.entities[1]);
    
    // Expect entity1 to be removed and only entity2 to exist
    assert(ecs.entities.find(1) != ecs.entities.end());
    assert(ecs.entities[1]->position.x == 0.0f); // Confirm removal was correct
}

void test_collision_resolution_accuracy() {
    ECS ecs(10, 10, 50.0f);

    ecs.insertEntity(1, vec2(100.0f, 0.0f), vec2(+2.0f, 0.0f), ShapeType::Circle, 0.01f, 1.0f);
    ecs.insertEntity(2, vec2(110.0f, 0.0f), vec2(-2.0f, 0.0f), ShapeType::Circle, 0.01f, 10.0f);

    ecs.tick_simulation_precise(5.0f);  // Simulate 5 seconds

    // ecs.printEntityVelocities();
    // ecs.printEntityPositions();

    Entity* entity1 = ecs.entities[1];
    Entity* entity2 = ecs.entities[2];

    // Expect velocities to be swapped based on mass ratio
    assert_apeq(entity1->velocity.x, -5.27273f);
    assert_apeq(entity2->velocity.x, -1.27273f);

    assert_apeq(entity1->position.x, 100-8.21803f);
    assert_apeq(entity2->position.x, 100+1.8218f);
}

void test_big_objects() {
    ECS ecs(100, 100, 10.0f);

    ecs.insertEntity(1, vec2(  0.0f   , 0.0f), vec2(+2.0f, 0.0f), ShapeType::Circle, 50.0f, 1.0f);
    ecs.insertEntity(2, vec2(110.0001f, 0.0f), vec2(-2.0f, 0.0f), ShapeType::Circle, 50.0f, 10.0f);

    ecs.printEntityFulls();
    ecs.tick_simulation_precise(5.0f);  // Simulate 5 seconds

    // ecs.printEntityVelocities();
    ecs.printEntityFulls();

    Entity* entity1 = ecs.entities[1];
    Entity* entity2 = ecs.entities[2];

    // Expect velocities to be swapped based on mass ratio
    assert_apeq(entity1->velocity.x, -5.27273f);
    assert_apeq(entity2->velocity.x, -1.27273f);
    
    assert_apeq(entity1->position.x, 000-8.18164f);
    assert_apeq(entity2->position.x, 100+1.8218f);
}

int main() {
    std::cout << "Running tests...\n";
    TEST(test_no_collision_movement());
    TEST(test_basic_collision());
    TEST(test_mass_impact_on_velocity());
    TEST(test_circle_circle_collision_time());
    TEST(test_entity_removal());
    TEST(test_collision_resolution_accuracy());
    TEST(test_big_objects());
    std::cout << "Tests completed.\n";
    
    const int screenWidth = 800;
    const int screenHeight = 450;
    InitWindow(screenWidth, screenHeight, "raylib [shapes] example - basic shapes drawing");
    SetTargetFPS(60);
    int cells = 128;
    ECS ecs(cells, cells, 800.0f / cells);
pl("\n")

    for (int i = 0; i < 20; i++) {
        float x = screenWidth/2.f;
        float y = (float(i)*screenHeight)/20;
        float radius = 5.0f;
        // float mass = (rand() % 50) / 10.0f + 1.0f;
        float mass = 10000.0f;
        ecs.insertEntity(i, vec2(x, y), vec2(0.01), ShapeType::Circle, radius, mass);
    }
    srand(2);
    int numEntities = 20000;
    for (int i = 0; i < numEntities; i++) {
        float x = rand() % screenWidth;
        float y = rand() % screenHeight;
        float vx = (rand() % 200 - 100) * 1.0f;
        float vy = (rand() % 200 - 100) * 1.0f;
        float radius = (rand() % 50) / 10.0f + 0.5f;
        // float radius = (rand() % 50) / 10.0f + 10.0f;
        // float mass = (rand() % 50) / 10.0f + 1.0f;
        float mass = radius*radius/100.f;
        ecs.insertEntity(i+20, vec2(x, y), vec2(vx, vy), ShapeType::Circle, radius, mass);
    }
        float x = 100;
        float y = 100;
        float vx = 50;
        float vy = 0;
        float radius = 50.0f;
        // float mass = (rand() % 50) / 10.0f + 1.0f;
        float mass = radius*radius/100.f;
        ecs.insertEntity(777, vec2(x, y), vec2(vx, vy), ShapeType::Circle, radius, mass);
        // x = 401;
        // vx = -50;
        // ecs.insertEntity(2, vec2(x, y), vec2(vx, vy), ShapeType::Circle, radius, mass);
        // x = 600;
        // vx = -80;
        // ecs.insertEntity(3, vec2(x, y), vec2(vx, vy), ShapeType::Circle, radius, mass);
        float total_energy = 0;
        for (auto& [id, entity] : ecs.entities) {
            total_energy += entity->mass * dot(entity->velocity,entity->velocity);
        }
        // std::cout << "******\n\n";
    // ecs.printEntityFulls();
    while (!WindowShouldClose()){
        ClearBackground(RAYWHITE);
        BeginDrawing();
            for (const auto& [id, entity] : ecs.entities) {
                DrawCircle(entity->position.x, entity->position.y, entity->radius, DARKBLUE);            
            }

            if(IsKeyDown(KEY_ENTER)){
            for (const auto& [id, entity] : ecs.entities) {
                float d = 100.1;
                if(entity->position.x < -10)            {entity->position.x = screenWidth/2.f;}
                if(entity->position.x > screenWidth+10) {entity->position.x = screenWidth/2.f;}
                if(entity->position.y < -10)             {entity->position.y = screenHeight/2.f;}
                if(entity->position.y > screenHeight+10) {entity->position.y = screenHeight/2.f;}

                if(entity->position.x < 0) {entity->velocity.x *= -1;}
                if(entity->position.x > screenWidth) {entity->velocity.x *= -1;}
                if(entity->position.y < 0) {entity->velocity.y *= -1;}
                if(entity->position.y > screenHeight) {entity->velocity.y *= -1;}
            }
                collisions.clear();
                if(IsKeyDown(KEY_LEFT_SHIFT)){
                    ecs.tick_simulation_precise(0.016);
                } else {
                    ecs.tick_simulation_fast(0.016);
                }
                // ecs.printEntityFulls();
            }
            // for (auto p : collisions){
            //     DrawCircle(p.x, p.y, 3, GREEN);   
            // }
            DrawFPS(10,10);
        EndDrawing();
    }
    // ecs.printEntityFulls();
    // ecs.printEntityProblems();
    // benchmark(ecs, 1);
        std::cout << "******\n\n";
    CloseWindow(); 

    std::cout << "total_energy is " << total_energy << "\n";
    for (auto& [id, entity] : ecs.entities) {
        total_energy -= entity->mass * dot(entity->velocity,entity->velocity);
    }
    std::cout << "energy difference is " << total_energy << "\n";
    std::cout << "collision counter is " << collisionCounter << "\n";

    return 0;
}