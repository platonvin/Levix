#include <chrono>
#include <iostream>
#include <glm/glm.hpp>
#include <unordered_map>
#include <vector>
#include <functional>
#include <cmath>

using namespace glm;

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

    int getCellIndex(vec2 position) {
        ivec2 xy = clamp(ivec2(position) / ivec2(cellSize), ivec2(0), size-1);
        return xy.x + xy.y * size.x;
    }

    void insertEntity(Entity* entity) {
        int index = getCellIndex(entity->position);
        entity->next = cells[index].head;
        cells[index].head = entity;
    }

    void removeEntity(Entity* entity) {
        int index = getCellIndex(entity->position);
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

    void moveEntity(Entity* entity, vec2 newPosition) {
        removeEntity(entity);
        entity->position = newPosition;
        insertEntity(entity);
    }

    void checkCollisions() {
        //for every cell
        for (int xx = 0; xx < size.x; xx++) {
        for (int yy = 0; yy < size.y; yy++) {
            ivec2 xxyy = ivec2(xx,yy);
            ivec2 low  = clamp(xxyy-1, ivec2(0), size-1);
            ivec2 high = clamp(xxyy+1, ivec2(0), size-1);

                int idx = xxyy.x + xxyy.y * size.x;
            Entity* entity = cells[idx].head;

            while (entity) {
                //for every nearby cell:
                for (int _x = low.x; _x <= high.x; _x++) {
                for (int _y = low.y; _y <= high.y; _y++) {
                        int idx = _x + _y * size.x;
                    Entity* other = cells[idx].head;
                    while (other) {
                        handleCollision(*entity, *other);
                        other = other->next;
                    }
                }}
                entity = entity->next;
            }
        }}
    }

    void handleCollision(Entity& entity1, Entity& entity2) {
        if (entity1.shape == ShapeType::Point && entity2.shape == ShapeType::Point) {
            return;  // No collision handling for two points
        }

        vec2 delta = entity2.position - entity1.position;
        float distance = length(delta);

        if (entity1.shape == ShapeType::Circle && entity2.shape == ShapeType::Circle) {
            // Circle vs Circle collision
            if (distance < entity1.radius + entity2.radius) {
                resolveCollision(entity1, entity2, delta, distance);
            }
            // std::cout << "Circle" << "\n";
        } else if ((entity1.shape == ShapeType::Point && entity2.shape == ShapeType::Circle) ||
                   (entity1.shape == ShapeType::Circle && entity2.shape == ShapeType::Point)) {
            // Point vs Circle collision
            float circleRadius = entity1.shape == ShapeType::Circle ? entity1.radius : entity2.radius;
            if (distance < circleRadius) {
                resolveCollision(entity1, entity2, delta, distance);
            }
            // std::cout << "Point" << "\n";
        } else if (entity1.shape == ShapeType::Square && entity2.shape == ShapeType::Square) {
            // AABB (Axis-Aligned Bounding Box) collision for squares
            if (abs(delta.x) < (entity1.width / 2 + entity2.width / 2) &&
                abs(delta.y) < (entity1.height / 2 + entity2.height / 2)) {
                resolveCollision(entity1, entity2, delta, distance);
            }
            // std::cout << "Square" << "\n";
        }
    }

private:
    // Resolve collision between two entities
    void resolveCollision(Entity& entity1, Entity& entity2, const vec2& delta, float distance) {
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

    void runCollisions() {
        spatialPartitioning.checkCollisions();
    }

    void _actual_tick_simulation(float deltaTime){
        for (auto& [id, entity] : entities) {
            entity->position += entity->velocity * deltaTime;
            spatialPartitioning.moveEntity(entity, entity->position);
        }
        runCollisions();
    }

    // Divided this way to make simulation more accurate
    void tick_simulation(float deltaTime) {
        float max_step = 0;
        for (auto& [id, entity] : entities) {
            max_step = max(max_step, length(entity->velocity * deltaTime));
        }
        if(max_step > spatialPartitioning.cellSize){
            int adjusted_step_count = ceil(max_step / spatialPartitioning.cellSize);        
            for (int i=0; i<adjusted_step_count; i++){
                _actual_tick_simulation(deltaTime / float(adjusted_step_count));
            }
            std::cout << adjusted_step_count << "\n";
        } else {
            _actual_tick_simulation(deltaTime);
        }
    }

    void printEntityPositions() {
        for (const auto& [id, entity] : entities) {
            std::cout << "Entity " << id << ": Position (" << entity->position.x << ", " << entity->position.y << ")\n";
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
    std::cout << duration.count() / numIterations * 1000<< " ms / iter\n";
}   

void test_no_collision_movement() {
    ECS ecs(10, 10, 50.0f);

    ecs.addEntity(1, vec2(0.0f, 0.0f), vec2(1.0f, 0.0f), ShapeType::Circle, 5.0f, 1.0f);
    ecs.tick_simulation(1.0f);  // Simulate 1 second of movement
    ecs.printEntityPositions();

    Entity* entity = ecs.entities[1];
    if (entity->position != vec2(1.0f, 0.0f)) {
        std::cout << "Test failed: Entity did not move as expected.\n";
    } else {
        std::cout << "Test passed: Entity moved correctly.\n";
    }
}

bool aprox_equal(float x, float y, float eps){
    return glm::abs(x-y) < abs(eps*x);
}
bool aprox_equal(float x, float y){
    return aprox_equal(x,y, 0.001);
}

void test_basic_collision() {
    ECS ecs(1, 1, 50.0f);

    ecs.addEntity(1, vec2( 0.0000f, 0.0f), vec2(+1.0f, 0.0f), ShapeType::Circle, 5.0f, .000001f);
    ecs.addEntity(2, vec2(10.0001f, 0.0f), vec2(-2.0f, 0.0f), ShapeType::Circle, 5.0f, 10000.0f);

    for (int x=0; x<240*25; x++)
        ecs.tick_simulation(0.001);

    // ecs.tick_simulation(0.1);
    
    ecs.printEntityPositions();

    Entity* entity1 = ecs.entities[1];
    Entity* entity2 = ecs.entities[2];
    std::cout << "entity1->velocity.x " << entity1->velocity.x << "\n";
    std::cout << "entity2->velocity.x " << entity2->velocity.x << "\n";
    if (entity1->velocity.x <= 0.0f && entity2->velocity.x >= 0.0f) {
        std::cout << "Test passed: Collision handled correctly.\n";
    } else {
        std::cout << "Test failed: Collision not handled as expected.\n";
    }
}

void test_mass_impact_on_velocity() {
    ECS ecs(10, 10, 50.0f);

    ecs.addEntity(1, vec2(0.0f, 0.0f), vec2(1.0f, 0.0f), ShapeType::Circle, 5.0f, 1.0f);
    ecs.addEntity(2, vec2(10.0f, 0.0f), vec2(-1.0f, 0.0f), ShapeType::Circle, 5.0f, 10.0f);

    ecs.tick_simulation(5.0f);  // Simulate 5 seconds, expect a collision halfway
    ecs.runCollisions();
    ecs.printEntityPositions();

    Entity* entity1 = ecs.entities[1];
    Entity* entity2 = ecs.entities[2];

    if (entity1->velocity.x < 0.0f && std::abs(entity2->velocity.x) < std::abs(entity1->velocity.x)) {
        std::cout << "Test passed: Mass correctly affected velocity.\n";
    } else {
        std::cout << "Test failed: Mass did not impact velocity as expected.\n";
    }
}

int main() {
    std::cout << "Running tests...\n";

    test_no_collision_movement();
    test_basic_collision();
    // test_mass_impact_on_velocity();

    std::cout << "Tests completed.\n";
    
    int sz = 100;
    ECS ecs(sz, sz, 5000.0f / sz);

    int numEntities = 5000;
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
    // benchmark(ecs, 1000);
        std::cout << "******\n\n";

    std::cout << "collisition counter is " << collisionCounter << "\n";

    return 0;
}

/*



*/