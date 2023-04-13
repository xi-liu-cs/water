using System.Collections;
using System.Collections.Generic;
using UnityEngine;



public class BoidManager : MonoBehaviour
{
    public class Boid {
        
        public float px = 0f;
        public float py = 0f;
        public float pz = 0f;
        public Vector3Int hashPosition;

        public float vx = 1f;
        public float vy = 1f;
        public float vz = 1f;
        
        public Vector3 position {
            get => new Vector3(this.px, this.py, this.pz);
            set {}
        }

        public Boid(Vector3 initPos) {
            this.px = initPos.x;
            this.py = initPos.y;
            this.pz = initPos.z;
            this.hashPosition = BoidManager.current.GetGridIndex(this.position);
        }
        public void Update() {
            // Get neighboring boids
            List<Boid> neighbors = new List<Boid>();
            List<Boid> closeNeighbors = new List<Boid>();
            BoidManager.current.GetNeighboringBoids(this ,out neighbors, out closeNeighbors);
            
            // local variables important for update
            float xpos_avg = 0f; 
            float ypos_avg = 0f;
            float zpos_avg = 0f;
            float xvel_avg = 0f;
            float yvel_avg = 0f;
            float zvel_avg = 0f;

            float close_dx = 0f;
            float close_dy = 0f;
            float close_dz = 0f;

            // First, we update `close_dx/dy/dz` based on boids within the protected range
            UpdateSeparation(closeNeighbors, ref close_dx, ref close_dy, ref close_dz);
            // Secondly, update `x/y/zpos_avg` and `x/y/zvel_avg` based on boids not in protected range but within visual range
            UpdateAlignment(neighbors, ref xvel_avg, ref yvel_avg, ref zvel_avg);
            UpdateCohesion(neighbors, ref xpos_avg, ref ypos_avg, ref zpos_avg);
            // Thirdly, update velocity based on if there were neighbor boids within visual range but not in protected range
            if (neighbors.Count > 0) {
                this.vx = this.vx 
                    + (xpos_avg - this.px) * BoidManager.current.centeringFactor 
                    + (xvel_avg - this.vx) * BoidManager.current.matchingFactor;
                this.vy = this.vy 
                    + (ypos_avg - this.py) * BoidManager.current.centeringFactor 
                    + (yvel_avg - this.vy) * BoidManager.current.matchingFactor;
                this.vz = this.vz 
                    + (zpos_avg - this.pz) * BoidManager.current.centeringFactor 
                    + (zvel_avg - this.vz) * BoidManager.current.matchingFactor;
            }
            // Fourthly, add avoidance contribution
            this.vx = this.vx + (close_dx * BoidManager.current.avoidFactor);
            this.vy = this.vy + (close_dy * BoidManager.current.avoidFactor);
            this.vz = this.vz + (close_dz * BoidManager.current.avoidFactor);
            // Fifthly, update based on edges and speed limits
            UpdateEdges();
            UpdateSpeedLimits();
            // Sixthly, update position
            this.px += this.vx;
            this.py += this.vy;
            this.pz += this.vz;
            // Lastly, update hash position
            Vector3Int newHashPosition = BoidManager.current.GetGridIndex(this.position);
            if (newHashPosition != this.hashPosition) {
                BoidManager.current.UpdateBoidGridPosition(this,this.hashPosition,newHashPosition);
            }
            this.hashPosition = newHashPosition;
            
        }
        public void UpdateSeparation(List<Boid> closeNeighbors, ref float close_dx, ref float close_dy, ref float close_dz) {
            if (closeNeighbors.Count == 0) return;
            foreach(Boid boid in closeNeighbors) {
                close_dx += this.position.x - boid.position.x;
                close_dy += this.position.y - boid.position.y;
                close_dz += this.position.z - boid.position.z;
            }
        }
        public void UpdateAlignment(List<Boid> neighbors, ref float xvel_avg, ref float yvel_avg, ref float zvel_avg) {
            if(neighbors.Count == 0) return;
            foreach(Boid boid in neighbors) {
                xvel_avg += boid.vx;
                yvel_avg += boid.vy;
                zvel_avg += boid.vz;
            }
            xvel_avg /= neighbors.Count;
            yvel_avg /= neighbors.Count;
            zvel_avg /= neighbors.Count;
        }
        public void UpdateCohesion(List<Boid> neighbors, ref float xpos_avg, ref float ypos_avg, ref float zpos_avg) {
            if (neighbors.Count == 0) return;
            foreach(Boid boid in neighbors) {
                xpos_avg += boid.px;
                ypos_avg += boid.py;
                zpos_avg += boid.pz;
            }
            xpos_avg /= neighbors.Count;
            ypos_avg /= neighbors.Count;
            zpos_avg /= neighbors.Count;
        }
        public void UpdateEdges() {
            if(this.px < BoidManager.current.minX) this.vx += BoidManager.current.turnFactor;
            if(this.px > BoidManager.current.maxX) this.vx -= BoidManager.current.turnFactor;
            if(this.py < BoidManager.current.minY) this.vy += BoidManager.current.turnFactor;
            if(this.py > BoidManager.current.maxY) this.vy -= BoidManager.current.turnFactor;
            if(this.pz < BoidManager.current.minZ) this.vz += BoidManager.current.turnFactor;
            if(this.pz > BoidManager.current.maxZ) this.vz -= BoidManager.current.turnFactor;
            return;
        }
        public void UpdateSpeedLimits() {
            float speed = Mathf.Sqrt( Mathf.Pow(this.vx,2) + Mathf.Pow(this.vy,2) + Mathf.Pow(this.vz,2) );
            if(speed > BoidManager.current.speedLimits.y) {
                this.vx = (this.vx / speed) * BoidManager.current.speedLimits.y;
                this.vy = (this.vy / speed) * BoidManager.current.speedLimits.y;
                this.vz = (this.vz / speed) * BoidManager.current.speedLimits.y;
            }
            if(speed < BoidManager.current.speedLimits.x) {
                this.vx = (this.vx / speed) * BoidManager.current.speedLimits.x;
                this.vy = (this.vy / speed) * BoidManager.current.speedLimits.x;
                this.vz = (this.vz / speed) * BoidManager.current.speedLimits.x;
            }
            return;
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
    
    [Range(1,20)]
    public int visualRange = 3;
    [Range(0,1)]
    public int protectedRange = 1;
    
    [Range(0f,0.5f)]
    public float turnFactor = 0.2f;
    [Range(0f,0.01f)]
    public float centeringFactor = 0.0005f;
    [Range(0f,0.1f)]
    public float avoidFactor = 0.05f;
    [Range(0f,0.1f)]
    public float matchingFactor = 0.05f;
    public Vector2 speedLimits = new Vector2(3f,6f);

    private Dictionary<Vector3Int,List<Boid>> grid = new Dictionary<Vector3Int,List<Boid>>();
    private List<Boid> boids = new List<Boid>();
    public float renderSize = 0.5f;

    public Color boidColor = Color.blue;
    public Color cellColor = Color.red;
    public bool drawBoids = true;
    public bool drawCurrentCells = true;

    void OnDrawGizmos() {
        
        Gizmos.color = Color.white;
        Gizmos.DrawWireCube(transform.position + dimensions/2f, dimensions);

        if(Application.isPlaying) {
            
            if(drawBoids) {
                Gizmos.color = boidColor;
                foreach(Boid boid in boids)
                    Gizmos.DrawSphere(boid.position,renderSize);
            }
            if(drawCurrentCells) {
                Color gridColor = cellColor;
                Vector3 divisionSize = new Vector3(
                    dimensions.x / numDivisions.x,
                    dimensions.y / numDivisions.y,
                    dimensions.z / numDivisions.z
                );
                foreach(KeyValuePair<Vector3Int,List<Boid>> kvp in grid) {
                    gridColor.a = Mathf.Clamp((float)kvp.Value.Count / ((float)numBoids/10f), 0f, 1f);
                    Gizmos.color = gridColor;
                    Vector3 divisionPos = new Vector3(
                        transform.position.x + kvp.Key.x*divisionSize.x,
                        transform.position.y + kvp.Key.y*divisionSize.y,
                        transform.position.z + kvp.Key.z*divisionSize.z
                    ) + divisionSize/2f;    
                    Gizmos.DrawCube(divisionPos,divisionSize);
                }
            }
            /*
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
            */
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
                                if(isClose) closeNeighbors.Add(neighborBoid);
                                else neighbors.Add(neighborBoid);
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
