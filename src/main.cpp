#include <chrono>
#include <iostream>
#include <glm/glm.hpp>
#include <unordered_map>
#include <vector>
#include <functional>
#include <cmath>

#define pl(x) std::cout << #x " " << x << "\n";

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

using namespace glm;

// Do you really need more than that?
enum class ShapeType { 
    Point, 
    Circle, 
    Square 
};

struct Entity {
    int id;
    vec2 position;
    vec2 velocity;
    ShapeType shape;
        float radius;
        float width;
        float height;
    float mass;
    Entity* next;
    std::function<void(Entity&, Entity&)> onCollision;
};

struct CollisionEvent {
    Entity* entity1;
    Entity* entity2;
    float timeUntilCollision;
};

struct SpatialCell {
    Entity* head;
};

static long int collisionCounter = 0;

class SpatialPartitioning {
public:
    // int gridWidth, gridHeight;
    ivec2 size; //cells_x; cells_y
    float cellSize;
    std::vector<SpatialCell> cells;

    SpatialPartitioning(int gridWidth, int gridHeight, float cellSize)
        : size(gridWidth, gridHeight), cellSize(cellSize),
          cells(gridWidth * gridHeight) {}
    ~SpatialPartitioning() {}

    // This function will return all cell indices that the entity overlaps with
    std::vector<int> getCellIndices(Entity* entity) {
        // Determine the bounding box
        ivec2 minCell;
        ivec2 maxCell;

        switch (entity->shape) {
            case (ShapeType::Circle): {
                vec2 minPos = entity->position - vec2(entity->radius, entity->radius);
                vec2 maxPos = entity->position + vec2(entity->radius, entity->radius);
                minCell = clamp(ivec2(minPos) / ivec2(cellSize), ivec2(0), size - 1);
                maxCell = clamp(ivec2(maxPos) / ivec2(cellSize), ivec2(0), size - 1);
                break;
            }
            case (ShapeType::Square): {
                vec2 minPos = entity->position - vec2(entity->width / 2, entity->height / 2);
                vec2 maxPos = entity->position + vec2(entity->width / 2, entity->height / 2);
                minCell = clamp(ivec2(minPos) / ivec2(cellSize), ivec2(0), size - 1);
                maxCell = clamp(ivec2(maxPos) / ivec2(cellSize), ivec2(0), size - 1);
                break;
            }
            case (ShapeType::Point):
                break;

        }

      // Get the min and max cell indices

      std::vector<int> cellIndices;
      // Collect all cells within the bounding box
      for (int x = minCell.x; x <= maxCell.x; ++x) {
        for (int y = minCell.y; y <= maxCell.y; ++y) {
          cellIndices.push_back(x + y * size.x);
        }
      }
      return cellIndices;
    }

    void insertEntity(Entity* entity) {
        std::vector<int> indices = getCellIndices(entity);
        for (int index : indices) {
            entity->next = cells[index].head;
            cells[index].head = entity;
        }
    }

    void removeEntity(Entity* entity) {
        std::vector<int> indices = getCellIndices(entity);
        for (int index : indices) {
            Entity** current = &cells[index].head;
            while (*current != nullptr) {
                if (*current == entity) {
                    *current = entity->next;
                    entity->next = nullptr;
                    break;
                }
                current = &(*current)->next;
            }
        }
    }

    void moveEntity(Entity* entity, vec2 newPosition) {
        removeEntity(entity);
        entity->position = newPosition;
        insertEntity(entity);
    }
};

struct ECS {
    std::unordered_map<int, Entity*> entities;
    SpatialPartitioning spatialPartitioning;

    ECS(int gridWidth, int gridHeight, float cellSize)
        : spatialPartitioning(gridWidth, gridHeight, cellSize) {}
    ~ECS() {}

    void addEntity(int id, vec2 position, vec2 velocity, ShapeType shape, float size, float mass) {
        Entity* entity = new Entity{id, position, velocity, shape, 0, 0, 0, mass, nullptr};
        if (shape == ShapeType::Circle) {
            entity->radius = size;
        } else if (shape == ShapeType::Square) {
            entity->width = size;
            entity->height = size;
        }
        entities[id] = entity;
        spatialPartitioning.insertEntity(entity);
    }

    void updateEntityPosition(int id, vec2 newPosition) {
        Entity* entity = entities[id];
        spatialPartitioning.moveEntity(entity, newPosition);
    }

    void processSimulation(float deltaTime) {
        float usedTime = 0.0f;

        // Continue finding and processing collisions until we've exhausted the delta time
        // this is like general physics engine approach but precise and painfully slow
        while (usedTime < deltaTime) {
            // pl("START ITER")
            CollisionEvent nextCollision;
            bool collides = findNextCollision(deltaTime - usedTime, &nextCollision);
            if(collides){
                // pl(nextCollision.timeUntilCollision)

                // Move all entities by the time until the next collision
                float collisionTime = nextCollision.timeUntilCollision;
                moveEntitiesByTime(collisionTime);

                // pl("AHTUNG COLLISION HAPPENED")
                // printEntityPositions();
                // pl("\n")

                // Handle the next collision
                if (nextCollision.entity1 && nextCollision.entity2) {
                    // handleCollision(*nextCollision.entity1, *nextCollision.entity2);
                        vec2 delta = nextCollision.entity2->position - nextCollision.entity1->position;
                    resolveCollision(*nextCollision.entity1, *nextCollision.entity2, delta);
                }
                usedTime += collisionTime;
            } {
                moveEntitiesByTime(deltaTime - usedTime);
                break;
            }

            // Accumulate the used time
            // pl("END ITER")
        }
    }

    void moveEntitiesByTime(float time) {
        for (auto& [id, entity] : entities) {
            entity->position += entity->velocity * time;
            spatialPartitioning.moveEntity(entity, entity->position);
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

    bool findNextCollision(float remainingTime, CollisionEvent* event) {
        *event = {0};
        float closest_collision_time = +INFINITY;
        bool any_collision_found = false;
        
        // Check every pair of nearby entities and calculate when they'll collide
        for (int xx = 0; xx < spatialPartitioning.size.x; xx++) {
            for (int yy = 0; yy < spatialPartitioning.size.y; yy++) {
                ivec2 xxyy = ivec2(xx, yy);
                ivec2 low  = clamp(xxyy - 1, ivec2(0), spatialPartitioning.size - 1);
                ivec2 high = clamp(xxyy + 1, ivec2(0), spatialPartitioning.size - 1);

                int idx = xxyy.x + xxyy.y * spatialPartitioning.size.x;
                Entity* entity = spatialPartitioning.cells[idx].head;

                while (entity) {
                    // Check every nearby cell
                    for (int _x = low.x; _x <= high.x; _x++) {
                        for (int _y = low.y; _y <= high.y; _y++) {
                            int idx = _x + _y * spatialPartitioning.size.x;
                            Entity* other = spatialPartitioning.cells[idx].head;

                            while (other) {
                                // Calculate time until collision between entity and other
                                //TODO move checks around
                                if(other != entity){ //reduntant check
                                    // pl(entity)
                                    // pl(other)
                                    float timeUntilCollision;
                                    bool collides = calculateCollisionTime(*entity, *other, &timeUntilCollision);
                                    // pl(collides)

                                    if(collides){
                                        any_collision_found = true;
                                        // pl(timeUntilCollision);
                                        if ((timeUntilCollision < closest_collision_time) and
                                            noneEqual<Entity*>({entity, other}, {event->entity1, event->entity2})) {
                                            closest_collision_time = timeUntilCollision;
                                            *event = {entity, other, timeUntilCollision};
                                            // pl(event->entity1)
                                            // pl(event->entity2)
                                        }
                                    }
                                }
                                other = other->next;
                            }
                        }
                    }
                    entity = entity->next;
                }
            }
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

    void tick_simulation(float deltaTime) {
        float max_step = 0;
        for (auto& [id, entity] : entities) {
            max_step = max(max_step, length(entity->velocity * deltaTime));
        }
        //TODO no branch
        // Divided this way to make simulation more accurate
        if(max_step > spatialPartitioning.cellSize){
            int adjusted_step_count = ceil(max_step / spatialPartitioning.cellSize);        
            for (int i=0; i<adjusted_step_count; i++){
                processSimulation(deltaTime / float(adjusted_step_count));
            }
            std::cout << adjusted_step_count << "\n";
        } else {
            processSimulation(deltaTime);
        }
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
};

void benchmark(ECS& ecs, int numIterations) {
    auto startTime = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < numIterations; ++i) {
        ecs.tick_simulation(0.016f);
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = endTime - startTime;

    std::cout << "Benchmark: " << numIterations << " iterations took "
              << duration.count() << " seconds.\n";
    std::cout << (duration.count() / numIterations) * 1000 << " ms / iter\n";
}   

void test_no_collision_movement() {
    ECS ecs(10, 10, 50.0f);

    ecs.addEntity(1, vec2(0.0f, 0.0f), vec2(1.0f, 0.0f), ShapeType::Circle, 5.0f, 1.0f);
    ecs.tick_simulation(1.0f);  // Simulate 1 second of movement
    // ecs.printEntityPositions();

    Entity* entity = ecs.entities[1];
    assert_apeq(entity->position.x, +1.0f);
}

void test_basic_collision() {
    ECS ecs(1, 1, +INFINITY);

    ecs.addEntity(1, vec2( 0.0000f, 0.0f), vec2(+1.0f, 0.0f), ShapeType::Circle, 1.0f, .001f);
    ecs.addEntity(2, vec2(10.0001f, 0.0f), vec2(-1.0f, 0.0f), ShapeType::Circle, 1.0f, 10000000.0f);

    ecs.tick_simulation(100000000.0);
    
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

    ecs.addEntity(1, vec2(0.0f, 0.0f), vec2(1.0f, 0.0f), ShapeType::Circle, 5.0f, 1.0f);
    ecs.addEntity(2, vec2(10.0f, 0.0f), vec2(-1.0f, 0.0f), ShapeType::Circle, 5.0f, 10.0f);

    ecs.tick_simulation(5.0f); // Simulate 5 seconds, expect a collision halfway
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

    ecs.addEntity(1, vec2(0.0f, 0.0f), vec2(2.0f, 0.0f), ShapeType::Circle, 5.0f, 1.0f);
    ecs.addEntity(2, vec2(20.0f, 0.0f), vec2(-2.0f, 0.0f), ShapeType::Circle, 5.0f, 1.0f);

    ecs.tick_simulation(2.5f);
    // Expect the entities to collide after 2.5 seconds
    Entity* entity1 = ecs.entities[1];
    Entity* entity2 = ecs.entities[2];
    assert_apeq(entity1->position.x, 5.0f);
    assert_apeq(entity2->position.x, 15.0f);

    // ecs.printEntityPositions();
}

void test_entity_removal() {
    ECS ecs(10, 10, 50.0f);

    ecs.addEntity(1, vec2(0.0f, 0.0f), vec2(1.0f, 0.0f), ShapeType::Circle, 5.0f, 1.0f);
    ecs.addEntity(2, vec2(10.0f, 0.0f), vec2(-1.0f, 0.0f), ShapeType::Circle, 5.0f, 1.0f);

    ecs.tick_simulation(0.0f);  // Simulate no movement

    ecs.spatialPartitioning.removeEntity(ecs.entities[1]);
    
    // Expect entity1 to be removed and only entity2 to exist
    assert(ecs.entities.find(1) != ecs.entities.end());
    assert(ecs.entities[1]->position.x == 0.0f); // Confirm removal was correct
}

void test_collision_resolution_accuracy() {
    ECS ecs(10, 10, 50.0f);

    ecs.addEntity(1, vec2(100.0f, 0.0f), vec2(+2.0f, 0.0f), ShapeType::Circle, 0.01f, 1.0f);
    ecs.addEntity(2, vec2(110.0f, 0.0f), vec2(-2.0f, 0.0f), ShapeType::Circle, 0.01f, 10.0f);

    ecs.tick_simulation(5.0f);  // Simulate 5 seconds

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

    ecs.addEntity(1, vec2(  0.0f   , 0.0f), vec2(+2.0f, 0.0f), ShapeType::Circle, 50.0f, 1.0f);
    ecs.addEntity(2, vec2(110.0001f, 0.0f), vec2(-2.0f, 0.0f), ShapeType::Circle, 50.0f, 10.0f);

    ecs.tick_simulation(5.0f);  // Simulate 5 seconds

    // ecs.printEntityVelocities();
    // ecs.printEntityPositions();

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
    test_no_collision_movement();
    test_basic_collision();
    test_mass_impact_on_velocity();
    test_circle_circle_collision_time();
    test_entity_removal();
    test_collision_resolution_accuracy();
    test_big_objects();
    std::cout << "Tests completed.\n";
    
    int cells = 500;
    ECS ecs(cells, cells, 6000.0f / cells);

    int numEntities = 1000;
    for (int i = 0; i < numEntities; ++i) {
        float x = rand() % 5000;
        float y = rand() % 5000;
        float vx = (rand() % 200 - 100) / 100.0f;
        float vy = (rand() % 200 - 100) / 100.0f;
        float radius = (rand() % 50) / 10.0f + 1.0f;
        float mass = (rand() % 50) / 10.0f + 1.0f;
        ecs.addEntity(i, vec2(x, y), vec2(vx, vy), ShapeType::Circle, radius, mass);
    }

        std::cout << "******\n\n";
    benchmark(ecs, 1);
        std::cout << "******\n\n";

    std::cout << "collision counter is " << collisionCounter << "\n";

    return 0;
}