using System.Collections;
using System.Collections.Generic;
using UnityEngine;



public class BoidManager : MonoBehaviour
{
    public class Boid {
        
        public Vector3 position;
        public Vector3Int hashPosition;
        public Vector3 velocity = Vector3.one;
        
        public Boid(Vector3 initPos) {
            this.position = initPos;
            this.hashPosition = BoidManager.current.GetGridIndex(this.position);
        }
        public void Update() {
            // Get neighboring boids
            List<Boid> neighbors = new List<Boid>();
            List<Boid> closeNeighbors = new List<Boid>();
            BoidManager.current.GetNeighboringBoids(this ,out neighbors, out closeNeighbors);
            
            // All update methods
            UpdateSeparation(closeNeighbors);
            UpdateAlignment(neighbors);
            UpdateCohesion(neighbors);
            UpdateEdges();
            UpdateSpeedLimits();

            // update position
            this.position += this.velocity;

            // Update hash position
            Vector3Int newHashPosition = BoidManager.current.GetGridIndex(this.position);
            if (newHashPosition != this.hashPosition) {
                BoidManager.current.UpdateBoidGridPosition(this,this.hashPosition,newHashPosition);
            }
            this.hashPosition = newHashPosition;
            
        }
        public void UpdateSeparation(List<Boid> closeNeighbors) {
            float close_dx = 0f, close_dy = 0f, close_dz = 0f;
            foreach(Boid boid in closeNeighbors) {
                close_dx += this.position.x - boid.position.x;
                close_dy += this.position.y - boid.position.y;
                close_dz += this.position.z - boid.position.z;
            }
            this.velocity = new Vector3(
                this.velocity.x + close_dx * BoidManager.current.avoidFactor,
                this.velocity.y + close_dy * BoidManager.current.avoidFactor,
                this.velocity.z + close_dz * BoidManager.current.avoidFactor
            );
            return;
        }
        public void UpdateAlignment(List<Boid> neighbors) {
            float xvel_avg = 0f, yvel_avg = 0f, zvel_avg = 0f;
            foreach(Boid boid in neighbors) {
                xvel_avg += boid.velocity.x;
                yvel_avg += boid.velocity.y;
                zvel_avg += boid.velocity.z;
            }
            if(neighbors.Count > 0) {
                xvel_avg /= neighbors.Count;
                yvel_avg /= neighbors.Count;
                zvel_avg /= neighbors.Count;
            }    
            this.velocity = new Vector3(
                this.velocity.x + (xvel_avg - this.velocity.x) * BoidManager.current.matchingFactor,
                this.velocity.y + (yvel_avg - this.velocity.y) * BoidManager.current.matchingFactor,
                this.velocity.z + (zvel_avg - this.velocity.z) * BoidManager.current.matchingFactor
            );
            return;
        }
        public void UpdateCohesion(List<Boid> neighbors) {
            float xpos_avg = 0f, ypos_avg = 0f, zpos_avg = 0f;
            foreach(Boid boid in neighbors) {
                xpos_avg += boid.position.x;
                ypos_avg += boid.position.y;
                zpos_avg += boid.position.z;
            }
            if (neighbors.Count > 0) {
                xpos_avg /= neighbors.Count;
                ypos_avg /= neighbors.Count;
                zpos_avg /= neighbors.Count;
            }
            this.velocity = new Vector3(
                this.velocity.x + (xpos_avg - this.position.x) * BoidManager.current.centeringFactor,
                this.velocity.y + (ypos_avg - this.position.y) * BoidManager.current.centeringFactor,
                this.velocity.z + (zpos_avg - this.position.z) * BoidManager.current.centeringFactor
            );
            return;
        }
        public void UpdateEdges() {
            float vx = this.velocity.x;
            float vy = this.velocity.y;
            float vz = this.velocity.z;
            if(this.position.x < BoidManager.current.minX) vx += BoidManager.current.turnFactor;
            if(this.position.x > BoidManager.current.maxX) vx -= BoidManager.current.turnFactor;
            if(this.position.y < BoidManager.current.minY) vy += BoidManager.current.turnFactor;
            if(this.position.y > BoidManager.current.maxY) vy -= BoidManager.current.turnFactor;
            if(this.position.z < BoidManager.current.minZ) vz += BoidManager.current.turnFactor;
            if(this.position.z > BoidManager.current.maxZ) vz -= BoidManager.current.turnFactor;
            this.velocity = new Vector3(vx,vy,vz);
        }
        public void UpdateSpeedLimits() {
            float speed = this.velocity.magnitude;
            if(speed > BoidManager.current.speedLimits.y) {
                this.velocity = (this.velocity / speed) * BoidManager.current.speedLimits.y;
            }
            if(speed < BoidManager.current.speedLimits.x) {
                this.velocity = (this.velocity / speed) * BoidManager.current.speedLimits.x;
            }
        }
    }

    public static BoidManager current;

    public int numBoids = 100;
    public Vector3 dimensions = new Vector3(200f,200f,200f);
    public Vector3Int numDivisions = new Vector3Int(10,10,10); 
    public float minX {
        get => transform.position.x;
        set {}
    }
    public float maxX {
        get => transform.position.x + dimensions.x;
        set {}
    }
    public float minY {
        get => transform.position.y;
        set {}
    }
    public float maxY {
        get => transform.position.y + dimensions.y;
        set {}
    }
    public float minZ {
        get => transform.position.z;
        set {}
    }
    public float maxZ {
        get => transform.position.z + dimensions.z;
        set {}
    }
    
    public int visualRange = 3;
    public int protectedRange = 1;
    
    public float turnFactor = 0.2f;
    public float centeringFactor = 0.0005f;
    public float avoidFactor = 0.05f;
    public float matchingFactor = 0.05f;
    public Vector2 speedLimits = new Vector2(3f,6f);
    
    public bool drawCurrentCells = true;
    public bool drawVisualCells = true;
    public bool drawCloseCells = true;

    private Dictionary<Vector3Int,List<Boid>> grid = new Dictionary<Vector3Int,List<Boid>>();
    private List<Boid> boids = new List<Boid>();
    public float renderSize = 0.5f;

    void OnDrawGizmos() {
        
        Gizmos.color = Color.white;
        Gizmos.DrawWireCube(transform.position + dimensions/2f, dimensions);

        if(Application.isPlaying) {
            List<Vector3Int> toPrintSections = new List<Vector3Int>();
            List<Vector3Int> visualSections = new List<Vector3Int>();
            List<Vector3Int> protectedSections = new List<Vector3Int>();
            Gizmos.color = Color.yellow;
            foreach(Boid boid in boids) {
                Gizmos.DrawSphere(boid.position,renderSize);
                Vector3Int hashIndexes = GetGridIndex(boid.position);
                if(!toPrintSections.Contains(hashIndexes)) toPrintSections.Add(hashIndexes);
                // Addding visual sections
                if (drawVisualCells) {
                    for(int x = hashIndexes.x-visualRange; x <= hashIndexes.x+visualRange; x++)
                        for(int y = hashIndexes.y-visualRange; y <= hashIndexes.y+visualRange;y++) 
                            for(int z = hashIndexes.z-visualRange; z <= hashIndexes.z+visualRange; z++) {
                                Vector3Int neighborIndex = new Vector3Int(x,y,z);
                                if(!visualSections.Contains(neighborIndex)) visualSections.Add(neighborIndex);
                            }
                }
                // Adding protected sections
                if (drawCurrentCells) {
                    for(int x = hashIndexes.x-protectedRange; x <= hashIndexes.x+protectedRange; x++)
                        for(int y = hashIndexes.y-protectedRange; y <= hashIndexes.y+protectedRange;y++) 
                            for(int z = hashIndexes.z-protectedRange; z <= hashIndexes.z+protectedRange; z++) {
                                Vector3Int neighborIndex = new Vector3Int(x,y,z);
                                if(!protectedSections.Contains(neighborIndex)) protectedSections.Add(neighborIndex);
                            }
                }
            }
            Vector3 divisionSize = new Vector3(
                dimensions.x / numDivisions.x,
                dimensions.y / numDivisions.y,
                dimensions.z / numDivisions.z
            );
            if (drawVisualCells) {
                Gizmos.color = Color.yellow;
                foreach(Vector3Int neighborIndex in visualSections) {
                    Vector3 divisionPos = new Vector3(
                        transform.position.x + neighborIndex.x*divisionSize.x,
                        transform.position.y + neighborIndex.y*divisionSize.y,
                        transform.position.z + neighborIndex.z*divisionSize.z
                    ) + divisionSize/2f;
                    Gizmos.DrawWireCube(divisionPos,divisionSize);
                }
            }
            if (drawCloseCells) {
                Gizmos.color = Color.red;
                foreach(Vector3Int neighborIndex in protectedSections) {
                    Vector3 divisionPos = new Vector3(
                        transform.position.x + neighborIndex.x*divisionSize.x,
                        transform.position.y + neighborIndex.y*divisionSize.y,
                        transform.position.z + neighborIndex.z*divisionSize.z
                    ) + divisionSize/2f;
                    Gizmos.DrawWireCube(divisionPos,divisionSize);
                }
            }
            if (drawCurrentCells) {
                Gizmos.color = Color.blue;
                foreach(Vector3Int actualIndex in toPrintSections) {
                    Vector3 divisionPos = new Vector3(
                        transform.position.x + actualIndex.x*divisionSize.x,
                        transform.position.y + actualIndex.y*divisionSize.y,
                        transform.position.z + actualIndex.z*divisionSize.z
                    ) + divisionSize/2f;
                    Gizmos.DrawWireCube(divisionPos,divisionSize);
                }
            }
            
        } else {
            if (!drawCurrentCells) return;
            Gizmos.color = Color.blue;
            for(int x = 0; x < numDivisions.x; x++)
                for(int y = 0; y < numDivisions.y; y++)
                    for(int z = 0; z < numDivisions.z; z++) {
                        Vector3 divisionSize = new Vector3(
                            dimensions.x / numDivisions.x,
                            dimensions.y / numDivisions.y,
                            dimensions.z / numDivisions.z
                        );
                        Vector3 divisionPos = new Vector3(
                            transform.position.x + x*divisionSize.x,
                            transform.position.y + y*divisionSize.y,
                            transform.position.z + z*divisionSize.z
                        );
                        Gizmos.DrawWireCube(divisionPos+divisionSize/2f,divisionSize);
                    }
        }
    }

    private void Awake() {
        current = this;
    }

    // Start is called before the first frame update
    private void Start() {
        // Initializing the grid
        InitializeGrid();

        // Initializing Boids
        InitializeBoids();

        foreach(Boid boid in boids) {
            if(!grid.ContainsKey(boid.hashPosition)) grid[boid.hashPosition] = new List<Boid>();
            grid[boid.hashPosition].Add(boid);
        }
    }

    // Update is called once per frame
    private void Update() {
        foreach(Boid boid in boids) {
            // Update boid manually
            boid.Update();
            // Update position of boid in grid
            
        }
         
    }

    private void InitializeGrid() {
        for(int x = 0; x < numDivisions.x; x++)
            for(int y = 0; y < numDivisions.y; y++)
                for(int z = 0; z < numDivisions.z; z++)
                    grid[new Vector3Int(x,y,z)] = new List<Boid>();
    }

    private void InitializeBoids() {
        for(int i = 0; i < numBoids; i++) {
            Vector3 initPosition = new Vector3(
                Random.Range(transform.position.x,transform.position.x+dimensions.x),
                Random.Range(transform.position.y,transform.position.y+dimensions.y),
                Random.Range(transform.position.z,transform.position.z+dimensions.z)
            );
            Boid boid = new Boid(initPosition);
            boids.Add(boid);
        }
    }

    private Vector3Int GetGridIndex(Vector3 pos) {
        int hashX = Mathf.FloorToInt(((pos.x-transform.position.x)/dimensions.x)*numDivisions.x);
        int hashY = Mathf.FloorToInt(((pos.y-transform.position.y)/dimensions.y)*numDivisions.y);
        int hashZ = Mathf.FloorToInt(((pos.z-transform.position.z)/dimensions.z)*numDivisions.z);
        return new Vector3Int(hashX, hashY, hashZ);
    }

    private void GetNeighboringBoids(Boid boid, out List<Boid> neighbors, out List<Boid> closeNeighbors) {
        neighbors = new List<Boid>();
        closeNeighbors = new List<Boid>();
        bool isClose = false;
        for(int x = boid.hashPosition.x - visualRange; x <= boid.hashPosition.x + visualRange; x++)
            for(int y = boid.hashPosition.y - visualRange; y <= boid.hashPosition.y + visualRange; y++)
                for(int z = boid.hashPosition.z - visualRange; z <= boid.hashPosition.z + visualRange; z++) {
                    Vector3Int hp = new Vector3Int(x,y,z);
                    if(!grid.ContainsKey(hp)) continue;
                    if(grid[hp].Count > 0) {
                        isClose = Mathf.Abs(x-boid.hashPosition.x)<=protectedRange && Mathf.Abs(y-boid.hashPosition.y)<=protectedRange && Mathf.Abs(z-boid.hashPosition.z)<=protectedRange;
                        foreach(Boid neighborBoid in grid[hp]) {
                            if(neighborBoid != boid) {
                                neighbors.Add(neighborBoid);
                                if(isClose) closeNeighbors.Add(neighborBoid);
                            }
                        }
                    }
                }
    }

    public void UpdateBoidGridPosition(Boid boid, Vector3Int oldHashPosition, Vector3Int newHashPosition) {
        if(grid[oldHashPosition].Contains(boid)) grid[oldHashPosition].Remove(boid);
        if(!grid.ContainsKey(newHashPosition)) grid[newHashPosition] = new List<Boid>();
        grid[newHashPosition].Add(boid);
    }
}
