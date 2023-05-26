using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

public class Grid3D : MonoBehaviour
{
    [Header("== WORLD PROPERTIES ==")]
    #if UNITY_EDITOR
    [Help("Control parameters concerning world space seting and general overhead.")]
    #endif    
        [SerializeField, Tooltip("Where in world space are we centering the simulation around?")]
        private Vector3 _origin = Vector3.zero;
        public Vector3 origin { get => _origin; set {} }
        public float[] originF { get => new float[3] { _origin.x, _origin.y, _origin.z }; set {} }

    [Header("== GRID PROPERTIES ==")]
    #if UNITY_EDITOR
    [Help("Control parameters concerning grid generation and dimensions.")]
    #endif
        [Tooltip("How big are our grid cells? Refers to length of the grid cell, assumes cubic cell.")] 
        public float gridCellSize = 2f;

        [SerializeField, Tooltip("How many cells should we expect per axis?")]
        private Vector3Int _numGridCellsPerAxis = new Vector3Int(10,10,10);
        public Vector3Int numGridCellsPerAxis { get => _numGridCellsPerAxis; set {} }
        public float[] numGridCellsPerAxisF { get => new float[3]{_numGridCellsPerAxis.x, _numGridCellsPerAxis.y, _numGridCellsPerAxis.z}; set {} }
        public int3 numGridCellsPerAxisI { get => new(_numGridCellsPerAxis.x, _numGridCellsPerAxis.y, _numGridCellsPerAxis.z); set {} }
        public int numGridCells { get => _numGridCellsPerAxis.x * _numGridCellsPerAxis.y * _numGridCellsPerAxis.z; set {} }

        [Tooltip("Automatically calculated world space length (per axis) of the simulation, based on `gridCellSize`. and `numGridCellsPerAxis`. ")]
        public Vector3 bounds { get => gridCellSize*(Vector3)_numGridCellsPerAxis; set {} }
        public float[] boundsF { get => new float[3]{bounds.x, bounds.y, bounds.z}; set {} }
        public Vector3Int numBufferCellsPerAxis = new Vector3Int(2,2,2);
        public Vector3 innerBounds { get => 
            gridCellSize * (Vector3)(_numGridCellsPerAxis - numBufferCellsPerAxis);
            set {}
        }
            

    [Header("== DEBUGGING ==")]
    public bool showBounds = false;
    public Color gridBoundsColor = Color.white;
    public bool showInnerBounds = false;
    public Color gridInnerBoundsColor = Color.yellow;
    public bool showGridAxes = false;
    public Color gridAxesColor = new Vector4(1f,1f,1f,0.5f);
    public bool debugGridStats = true;

    void OnDrawGizmos() {
        DrawGridAxes();
        DrawGridBounds();
        DrawGridInnerBounds();
    }

    public void DrawGridBounds(bool forceShow = false) {
        if (!forceShow && !showBounds) return;
        Gizmos.color = gridBoundsColor;
        Gizmos.DrawWireCube(_origin, bounds);
    }

    public void DrawGridInnerBounds(bool forceShow = false) {
        if (!forceShow && !showInnerBounds) return;
        Gizmos.color = gridInnerBoundsColor;
        Gizmos.DrawWireCube(_origin, innerBounds);
    }

    public void DrawGridAxes(bool forceShow = false) {
        if (!forceShow && !showGridAxes) return;
        Gizmos.color = gridAxesColor;
        Vector3 p1, p2, p3;
        for(int z = 0; z <= _numGridCellsPerAxis.z; z++) {
            p1 = Grid3DHelpers.GetGridCellWorldPositionFromXYZIndices(numGridCellsPerAxis, _origin, gridCellSize, 0, 0, z, false);
            p2 = Grid3DHelpers.GetGridCellWorldPositionFromXYZIndices(numGridCellsPerAxis, _origin, gridCellSize, numGridCellsPerAxis.x, 0, z, false);
            p3 = Grid3DHelpers.GetGridCellWorldPositionFromXYZIndices(numGridCellsPerAxis, _origin, gridCellSize, 0, numGridCellsPerAxis.y, z, false);
            Gizmos.DrawLine(p1,p2);
            Gizmos.DrawLine(p1,p3);
        }
        for(int x = 0; x <= _numGridCellsPerAxis.x; x++) {
            p1 = Grid3DHelpers.GetGridCellWorldPositionFromXYZIndices(numGridCellsPerAxis, _origin, gridCellSize, x, 0, 0, false);
            p2 = Grid3DHelpers.GetGridCellWorldPositionFromXYZIndices(numGridCellsPerAxis, _origin, gridCellSize, x, 0, numGridCellsPerAxis.z, false);
            Gizmos.DrawLine(p1,p2);
        }
        for(int y = 0; y <= _numGridCellsPerAxis.y; y++) {
            p1 = Grid3DHelpers.GetGridCellWorldPositionFromXYZIndices(numGridCellsPerAxis, _origin, gridCellSize, 0, y, 0, false);
            p2 = Grid3DHelpers.GetGridCellWorldPositionFromXYZIndices(numGridCellsPerAxis, _origin, gridCellSize, 0, y, numGridCellsPerAxis.z, false );
            Gizmos.DrawLine(p1,p2);
        }
    }

    void Awake() {
        if (!debugGridStats) return;
        Debug.Log($"GRID: [{_numGridCellsPerAxis.ToString()}]");
    }
}


public class Grid3DHelpers : MonoBehaviour
{
    // Get the projected index of a grid cell based on the xyz indices, given dimensions (# of cells along each axis)
    // The smallest indice range is the X. 
    // Moving along the Y axis within the same Z index, we can move by adding/subtracting X cells
    // Moving along the Z axis, we move by adding/subtracting X cells * Y cells
    public static int GetProjectedGridIndexFromXYZ(Vector3Int numGridCellsPerAxis, Vector3Int xyz) {
        return xyz.x + (numGridCellsPerAxis.x * xyz.y) + (numGridCellsPerAxis.x * numGridCellsPerAxis.y * xyz.z);
    }
    public static int GetProjectedGridIndexFromXYZ(Vector3Int numGridCellsPerAxis, int3 xyz) {
        return xyz[0] + (numGridCellsPerAxis.x * xyz[1]) + (numGridCellsPerAxis.x * numGridCellsPerAxis.y * xyz[2]);
    }
    public static int GetProjectedGridIndexFromXYZ(Vector3Int numGridCellsPerAxis, int x, int y, int z) {
        return x + (numGridCellsPerAxis.x * y) + (numGridCellsPerAxis.x * numGridCellsPerAxis.y * z);
    }

    // Given a boid, get the projected grid index from the boid's current world position, given bounds, global grid cell size, and dimensions
    public static int GetProjectedGridIndexFromGivenPosition(Vector3Int numGridCellsPerAxis, Vector3 origin, float gridCellSize, Vector3 position) {
        // Convert the given position to XYZ
        return GetProjectedGridIndexFromXYZ(numGridCellsPerAxis, GetGridXYZIndices(numGridCellsPerAxis, origin, gridCellSize, position));
    }

    // Get the XYZ Indices of a world position, given the bounds and the global size of a grid cell. 
    // The minimum bound is expected to be -bounds._ / 2f
    // Agnostic to negative/positive bound positions in world space. Will succeed no matter what.
    public static Vector3Int GetGridXYZIndices(Vector3Int numGridCellsPerAxis, Vector3 origin, float gridCellSize, Vector3 position) {
        return new Vector3Int(
            Mathf.FloorToInt((position.x - (origin.x - (numGridCellsPerAxis.x*gridCellSize)/2f))/gridCellSize),
            Mathf.FloorToInt((position.y - (origin.y - (numGridCellsPerAxis.y*gridCellSize)/2f))/gridCellSize),
            Mathf.FloorToInt((position.z - (origin.z - (numGridCellsPerAxis.z*gridCellSize)/2f))/gridCellSize)
        );
    }
    // Get the XYZ Indices of a world position, given the bounds and the global size of a grid cell. 
    // The minimum bound is expected to be -bounds._ / 2f
    public static int[] GetInt3GridXYZIndices(Vector3Int numGridCellsPerAxis, Vector3 origin, float gridCellSize, Vector3 position) {
        Vector3Int xyz = GetGridXYZIndices(numGridCellsPerAxis, origin, gridCellSize, position);
        return new int[3]{ xyz.x, xyz.y, xyz.z };
    }

    // Get the world position of a grid cell, given bounds and a global cell size
    // The minimum bound is expected to be -bounds._/2f
    public static Vector3 GetGridCellWorldPositionFromXYZIndices(Vector3Int numGridCellsPerAxis, Vector3 origin, float gridCellSize, Vector3Int xyz, bool getCenter = true) {
        float c = (getCenter) ? (gridCellSize/2f) : 0f;
        return new Vector3(
            (origin.x - ((numGridCellsPerAxis.x*gridCellSize)/2f)) + (xyz.x * gridCellSize) + c,
            (origin.y - ((numGridCellsPerAxis.y*gridCellSize)/2f)) + (xyz.y * gridCellSize) + c,
            (origin.z - ((numGridCellsPerAxis.z*gridCellSize)/2f)) + (xyz.z * gridCellSize) + c
        );
    }
    public static Vector3 GetGridCellWorldPositionFromXYZIndices(Vector3Int numGridCellsPerAxis, Vector3 origin, float gridCellSize, int x, int y, int z, bool getCenter = true) {
        float c = (getCenter) ? (gridCellSize/2f) : 0f;
        return new Vector3(
            (origin.x - ((numGridCellsPerAxis.x*gridCellSize)/2f)) + (x * gridCellSize) + c,
            (origin.y - ((numGridCellsPerAxis.y*gridCellSize)/2f)) + (y * gridCellSize) + c,
            (origin.z - ((numGridCellsPerAxis.z*gridCellSize)/2f)) + (z * gridCellSize) + c
        );
    }

    // Get the world position of a grid cell based on another (ex. boid's) world position, given bounds an a global grid cell size
    public static Vector3 GetGridCellWorldPositionFromGivenPosition(Vector3Int numGridCellsPerAxis, Vector3 origin, float gridCellSize, Vector3 position, bool getCenter = true) {
        Vector3Int xyz = GetGridXYZIndices(numGridCellsPerAxis, origin, gridCellSize, position);
        return GetGridCellWorldPositionFromXYZIndices(numGridCellsPerAxis, origin, gridCellSize, xyz, getCenter);
    }

}

