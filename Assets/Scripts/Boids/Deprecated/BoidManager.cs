using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Random = UnityEngine.Random;
using System.Linq;

public class BoidManager : MonoBehaviour
{

    public ComputeShader boidManagerShader;
    public ComputeBuffer boidsBuffer;
    public int boidUpdateKernel;
    public struct Boid2 {
        public float px;
        public float py;
        public float pz;
        public Vector3Int hashPosition;

        public float vx;
        public float vy;
        public float vz;
    };
    Boid2[] boids2;

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
            List<Vector3> farObstacles = new List<Vector3>();
            List<Vector3> closeObstacles = new List<Vector3>();
            BoidManager.current.GetNeighboringBoidsAndObstacles(
                this,
                out neighbors, 
                out closeNeighbors,
                out farObstacles,
                out closeObstacles
            );
            
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
            float far_obstacle_dx = 0f;
            float far_obstacle_dy = 0f;
            float far_obstacle_dz = 0f;
            float close_obstacle_dx = 0f;
            float close_obstacle_dy = 0f;
            float close_obstacle_dz = 0f;

            // First, we update `close_dx/dy/dz` based on boids within the protected range
            UpdateSeparation(closeNeighbors, ref close_dx, ref close_dy, ref close_dz);
            // Secondly, we update `far_obstacle_dx/y/z` and `close_obstacle_dx/y/z` based on detected obstacles
            UpdateObstacles(farObstacles, ref far_obstacle_dx, ref far_obstacle_dy, ref far_obstacle_dz);
            UpdateObstacles(closeObstacles, ref close_obstacle_dx, ref close_obstacle_dy, ref close_obstacle_dz);
            // Thirdly, update `x/y/zpos_avg` and `x/y/zvel_avg` based on boids not in protected range but within visual range
            UpdateAlignment(neighbors, ref xvel_avg, ref yvel_avg, ref zvel_avg);
            UpdateCohesion(neighbors, ref xpos_avg, ref ypos_avg, ref zpos_avg);
            // Fourthly, update velocity based on if there were neighbor boids within visual range but not in protected range
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
            // Fifthly, add avoidance contribution
            this.vx = this.vx 
                + (close_dx * BoidManager.current.avoidFactor) 
                + (far_obstacle_dx * BoidManager.current.obstacleFactor * 0.5f)
                + (close_obstacle_dx * BoidManager.current.obstacleFactor);
            this.vy = this.vy 
                + (close_dy * BoidManager.current.avoidFactor)
                + (far_obstacle_dy * BoidManager.current.obstacleFactor * 0.5f)
                + (close_obstacle_dy * BoidManager.current.obstacleFactor);
            this.vz = this.vz 
                + (close_dz * BoidManager.current.avoidFactor)
                + (far_obstacle_dz * BoidManager.current.obstacleFactor * 0.5f)
                + (close_obstacle_dz * BoidManager.current.obstacleFactor);
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
        public void UpdateObstacles(List<Vector3> obstacles, ref float obstacle_dx, ref float obstacle_dy, ref float obstacle_dz) {
            if(obstacles.Count == 0) return;
            foreach(Vector3 norm in obstacles) {
                obstacle_dx += norm.x;
                obstacle_dy += norm.y;
                obstacle_dz += norm.z;
            }
            obstacle_dx /= obstacles.Count;
            obstacle_dy /= obstacles.Count;
            obstacle_dz /= obstacles.Count;
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
        }
    }

    [Serializable]
    public class Obstacle {
        public GameObject gameObject;
        public bool isStatic = true;
        private Transform transform;
        private MeshFilter filter;
        public Vector3 prevPosition;
        public Quaternion prevRotation;
        public Vector3 prevScale;
        private Dictionary<Vector3Int, List<Vector3>> normalsDict = new Dictionary<Vector3Int, List<Vector3>>();

        public void Initialize() {
            this.prevPosition = this.gameObject.transform.position;
            this.transform = this.gameObject.transform;
            this.filter = this.gameObject.GetComponent<MeshFilter>();
            this.CalculateNormals();
        }
        public void Update() {
            if(filter == null || isStatic) return;
            bool changed = (
                this.prevPosition != this.transform.position
                || this.prevRotation != this.transform.rotation
                || this.prevScale != this.transform.localScale
            );
            if(changed) {
                // Our obstacle moved. That means that we need to recalculate the normals
                // First, let's update BoidManager so that the normals associated with each Vector3Int hashPosition is removed
               this.CalculateNormals();
            }
            this.prevPosition = this.transform.position;
            this.prevRotation = this.transform.rotation;
            this.prevScale = this.transform.localScale;
        }
        public void CalculateNormals() {
            Vector3[] verts = this.filter.mesh.vertices;
		    Vector4[] tangents = this.filter.mesh.tangents;
		    Vector3[] norms = this.filter.mesh.normals;
            Dictionary<Vector3Int, List<Vector3>> newNormalsDict = new Dictionary<Vector3Int, List<Vector3>>();
            for(int i = 0; i < verts.Length; i++) {
                Vector3 p = this.transform.TransformPoint(verts[i]);
                Vector3Int hp = BoidManager.current.GetGridIndex(p);
                Vector3 normal = this.transform.rotation * norms[i];
                if(!newNormalsDict.ContainsKey(hp)) newNormalsDict[hp] = new List<Vector3>();
                newNormalsDict[hp].Add(normal);
            }
            BoidManager.current.UpdateObstaclesGridVector(this.normalsDict,newNormalsDict);
            this.normalsDict = newNormalsDict;
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
    public float obstacleFactor = 0.05f;
    [Range(0f,0.1f)]
    public float matchingFactor = 0.05f;
    public Vector2 speedLimits = new Vector2(3f,6f);

    private Dictionary<Vector3Int,List<Boid>> grid = new Dictionary<Vector3Int,List<Boid>>();
    private Dictionary<Vector3Int,List<Vector3>> obstacleGrid = new Dictionary<Vector3Int, List<Vector3>>();
    private List<Boid> boids = new List<Boid>();
    public float renderSize = 0.5f;

    public Color boidColor = Color.blue;
    public Color cellColor = Color.red;
    public bool drawBoids = true;
    public bool drawCurrentCells = true;
    public bool drawObstacles = true;

    public List<Obstacle> obstacles = new List<Obstacle>();

    void OnDrawGizmos() {
        
        Gizmos.color = Color.white;
        Gizmos.DrawWireCube(transform.position + dimensions/2f, dimensions);

        if(Application.isPlaying) {
            
            if(drawBoids) {
                Gizmos.color = boidColor;
                foreach(Boid boid in boids)
                    Gizmos.DrawSphere(boid.position,renderSize);
            }
            Vector3 divisionSize = new Vector3(
                dimensions.x / numDivisions.x,
                dimensions.y / numDivisions.y,
                dimensions.z / numDivisions.z
            );
            if(drawCurrentCells) {
                Color gridColor = cellColor;
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

            if(drawObstacles && obstacles.Count > 0) {
                Gizmos.color = Color.white;
                foreach(KeyValuePair<Vector3Int,List<Vector3>> kvp in obstacleGrid) {
                    if (kvp.Value.Count == 0) continue;
                    Vector3 cellPos = new Vector3(
                        transform.position.x + kvp.Key.x*divisionSize.x,
                        transform.position.y + kvp.Key.y*divisionSize.y,
                        transform.position.z + kvp.Key.z*divisionSize.z
                    ) + divisionSize/2f;
                    Vector3 norm = Vector3.zero;
                    foreach(Vector3 normSmall in kvp.Value) {
                        norm += normSmall;
                    }
                    Gizmos.DrawLine(cellPos, cellPos + norm.normalized);
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

        InitializeGrid();
        InitializeBoids2();
        InitializeValues();
        InitializeBuffers();
    }

    private void InitializeBoids2() {
        boids2 = new Boid2[numBoids];
        for(int i = 0; i < numBoids; i++) {
            Vector3 initPosition = new Vector3(
                Random.Range(transform.position.x,transform.position.x+dimensions.x),
                Random.Range(transform.position.y,transform.position.y+dimensions.y),
                Random.Range(transform.position.z,transform.position.z+dimensions.z)
            );
            Vector3 initVel = new Vector3(
                Random.Range(-1f,1f),
                Random.Range(-1f,1f),
                Random.Range(-1f,1f)
            ).normalized;
            boids2[i] = new Boid2 {
                px = initPosition.x,
                py = initPosition.y,
                pz = initPosition.z,
                hashPosition = GetGridIndex(initPosition),
                vx = initVel.x,
                vy = initVel.y,
                vz = initVel.z
            };
        }
    }

    private void InitializeValues() {
        boidManagerShader.SetInt("numParticles",numBoids);
    }

    private void InitializeKernels() {
        boidUpdateKernel = boidManagerShader.FindKernel("UpdateBoids");
    }

    private void InitializeBuffers() {
        boidsBuffer = new ComputeBuffer(numBoids, sizeof(float)*6 + sizeof(int)*3);
        boidsBuffer.SetData(boids2);

        boidManagerShader.SetBuffer(boidUpdateKernel, "boids", boidsBuffer);
    }

    // Start is called before the first frame update
    private void Start() {
        // Initializing the grid
        InitializeGrid();

        // Initializing Boids
        InitializeBoids();

        // Initializing Obstacles
        InitializeObstacles();

        foreach(Boid boid in boids) {
            if(!grid.ContainsKey(boid.hashPosition)) grid[boid.hashPosition] = new List<Boid>();
            grid[boid.hashPosition].Add(boid);
        }
    }

    // Update is called once per frame
    /*
    private void Update() {
        foreach(Boid boid in boids) {
            // Update boid manually
            boid.Update();            
        }
        foreach(Obstacle obs in obstacles) {
            obs.Update();
        }
    }
    */

    private void Update() {
        // We first update our grid - hashing the positions of the boids, then producing a new array where boid ID is key and neighbors is 
    }

    private void InitializeGrid() {
        for(int x = 0; x < numDivisions.x; x++)
            for(int y = 0; y < numDivisions.y; y++)
                for(int z = 0; z < numDivisions.z; z++) {
                    Vector3Int p = new Vector3Int(x,y,z);
                    grid[p] = new List<Boid>();
                    //obstacleGrid[p] = new List<Vector3>();
                }
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

    private void InitializeObstacles() {
        foreach(Obstacle obstacle in obstacles) {
            obstacle.Initialize();
        }
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
    private void GetNeighboringBoidsAndObstacles(
        Boid boid, 
        out List<Boid> neighbors, 
        out List<Boid> closeNeighbors,
        out List<Vector3> farObstacles,
        out List<Vector3> nearbyObstacles
    ) {
        neighbors = new List<Boid>();
        closeNeighbors = new List<Boid>();
        farObstacles = new List<Vector3>();
        nearbyObstacles = new List<Vector3>();
        bool isClose;

        for(int x = boid.hashPosition.x - visualRange; x <= boid.hashPosition.x + visualRange; x++)
            for(int y = boid.hashPosition.y - visualRange; y <= boid.hashPosition.y + visualRange; y++)
                for(int z = boid.hashPosition.z - visualRange; z <= boid.hashPosition.z + visualRange; z++) {
                    Vector3Int hp = new Vector3Int(x,y,z);
                    isClose = Mathf.Abs(x-boid.hashPosition.x)<=protectedRange && Mathf.Abs(y-boid.hashPosition.y)<=protectedRange && Mathf.Abs(z-boid.hashPosition.z)<=protectedRange;
                    if(grid.ContainsKey(hp) && grid[hp].Count > 0) {
                        foreach(Boid neighborBoid in grid[hp]) {
                            if(neighborBoid == boid) continue;
                            if(isClose) closeNeighbors.Add(neighborBoid);
                            else neighbors.Add(neighborBoid);
                        }
                    }
                    if(obstacleGrid.ContainsKey(hp) && obstacleGrid[hp].Count > 0) {
                        foreach(Vector3 norm in obstacleGrid[hp]) {
                            if(isClose) nearbyObstacles.Add(norm);
                            else farObstacles.Add(norm);
                        }
                    }
                }
    }

    public void UpdateBoidGridPosition(Boid boid, Vector3Int oldHashPosition, Vector3Int newHashPosition) {
        if(grid[oldHashPosition].Contains(boid)) grid[oldHashPosition].Remove(boid);
        if(!grid.ContainsKey(newHashPosition)) grid[newHashPosition] = new List<Boid>();
        grid[newHashPosition].Add(boid);
    }
    public void UpdateObstaclesGridVector(Dictionary<Vector3Int, List<Vector3>> oldDict, Dictionary<Vector3Int, List<Vector3>> newDict) {
        foreach(KeyValuePair<Vector3Int, List<Vector3>> kvp in oldDict) {
            if(obstacleGrid.ContainsKey(kvp.Key) && obstacleGrid[kvp.Key].Count > 0) {
                foreach(Vector3 norm in kvp.Value) {
                    obstacleGrid[kvp.Key].Remove(norm);
                }
            }
        }
        foreach(KeyValuePair<Vector3Int, List<Vector3>> kvp in newDict) {
            if(!obstacleGrid.ContainsKey(kvp.Key)) obstacleGrid[kvp.Key] = new List<Vector3>();
            foreach(Vector3 norm in kvp.Value) {
                obstacleGrid[kvp.Key].Add(norm);
            }
        }
    }
}
