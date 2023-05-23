using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

public class Boid3DHelpers : MonoBehaviour
{
    // Get the projected index of a grid cell based on the xyz indices, given dimensions (# of cells along each axis)
    // The smallest indice range is the X. 
    // Moving along the Y axis within the same Z index, we can move by adding/subtracting X cells
    // Moving along the Z axis, we move by adding/subtracting X cells * Y cells
    public static int GetProjectedGridIndexFromXYZ(Vector3Int dimensions, Vector3Int xyz) {
        return (dimensions.x * dimensions.y * xyz.z) + (dimensions.x * xyz.y) + xyz.x;
    }
    public static int GetProjectedGridIndexFromXYZ(Vector3Int dimensions, int3 xyz) {
        return (dimensions.x * dimensions.y * xyz[2]) + (dimensions.x * xyz[1]) + xyz[0];
    }
    public static int GetProjectedGridIndexFromXYZ(Vector3Int dimensions, int x, int y, int z) {
        return (dimensions.x * dimensions.y * z) + (dimensions.x * y) + x;
    }

    // Given a boid, get the projected grid index from the boid's current world position, given bounds, global grid cell size, and dimensions
    public static int GetProjectedGridIndexFromGivenPosition(Vector3Int dimensions, Vector3 origin, float gridCellSize, Vector3 position) {
        // Convert the given position to XYZ
        return GetProjectedGridIndexFromXYZ(dimensions, GetGridXYZIndices(dimensions, origin, gridCellSize, position));
    }
    // Given a boid, get the projected grid index from the boid's current world position, given bounds, variable grid cell size, and dimensions
    public static int GetProjectedGridIndexFromGivenPosition(Vector3Int dimensions, Vector3 origin, Vector3 gridCellSizes, Vector3 position) {
        // Convert the given position to XYZ
        return GetProjectedGridIndexFromXYZ(dimensions, GetGridXYZIndices(dimensions, origin, gridCellSizes, position));
    }

    // Get the XYZ Indices of a world position, given the bounds and the global size of a grid cell. 
    // The minimum bound is expected to be -bounds._ / 2f
    public static Vector3Int GetGridXYZIndices(Vector3Int dimensions, Vector3 origin, float gridCellSize, Vector3 position) {
        // min bounds = 
        return new Vector3Int(
            Mathf.FloorToInt((position.x - (origin.x - (dimensions.x*gridCellSize)/2f))/gridCellSize),
            Mathf.FloorToInt((position.y - (origin.y - (dimensions.y*gridCellSize)/2f))/gridCellSize),
            Mathf.FloorToInt((position.z - (origin.z - (dimensions.z*gridCellSize)/2f))/gridCellSize)
        );
    }
    // Get the XYZ Indices of a world position, given the bounds and the global size of a grid cell. 
    // The minimum bound is expected to be -bounds._ / 2f
    public static int3 GetInt3GridXYZIndices(Vector3Int dimensions, Vector3 origin, float gridCellSize, Vector3 position) {
        // min bounds = 
        return new(
            Mathf.FloorToInt((position.x - (origin.x - (dimensions.x*gridCellSize)/2f))/gridCellSize),
            Mathf.FloorToInt((position.y - (origin.y - (dimensions.y*gridCellSize)/2f))/gridCellSize),
            Mathf.FloorToInt((position.z - (origin.z - (dimensions.z*gridCellSize)/2f))/gridCellSize)
        );
    }

    // Get the XYZ Indices of a world position, given the bounds and the variable size of a grid cell. 
    // The minimum bound is expected to be -bounds._ / 2f
    public static Vector3Int GetGridXYZIndices(Vector3Int dimensions, Vector3 origin, Vector3 gridCellSizes, Vector3 position) {
        return new Vector3Int(
            Mathf.FloorToInt((position.x - (origin.x - (dimensions.x*gridCellSizes.x)/2f))/gridCellSizes.x),
            Mathf.FloorToInt((position.y - (origin.y - (dimensions.y*gridCellSizes.y)/2f))/gridCellSizes.y),
            Mathf.FloorToInt((position.z - (origin.z - (dimensions.z*gridCellSizes.z)/2f))/gridCellSizes.z)
        );
        /*
        return new Vector3Int(
            Mathf.FloorToInt(position.x - (-bounds.x/2f)/gridCellSizes.x),
            Mathf.FloorToInt(position.y - (-bounds.y/2f)/gridCellSizes.y),
            Mathf.FloorToInt(position.z - (-bounds.z/2f)/gridCellSizes.z)
        );
        */
    }

    // Get the world position of a grid cell, given bounds and a global cell size
    // The minimum bound is expected to be -bounds._/2f
    public static Vector3 GetGridCellWorldPositionFromXYZIndices(Vector3Int dimensions, Vector3 origin, float gridCellSize, Vector3Int xyz) {
        return new Vector3(
            (origin.x - ((dimensions.x*gridCellSize)/2f)) + (xyz.x * gridCellSize) + (gridCellSize/2f),
            (origin.y - ((dimensions.y*gridCellSize)/2f)) + (xyz.y * gridCellSize) + (gridCellSize/2f),
            (origin.z - ((dimensions.z*gridCellSize)/2f)) + (xyz.z * gridCellSize) + (gridCellSize/2f)
        );
    }

    // Get the world position of a grid cell, given bounds and a variable cell size
    // The minimum bound is expected to be -bounds._/2f
    public static Vector3 GetGridCellWorldPositionFromXYZIndices(Vector3Int dimensions, Vector3 origin, Vector3 gridCellSizes, Vector3Int xyz) {
        return new Vector3(
            (origin.x - ((dimensions.x*gridCellSizes.x)/2f)) + (xyz.x * gridCellSizes.x) + (gridCellSizes.x/2f),
            (origin.y - ((dimensions.y*gridCellSizes.y)/2f)) + (xyz.y * gridCellSizes.y) + (gridCellSizes.y/2f),
            (origin.z - ((dimensions.z*gridCellSizes.z)/2f)) + (xyz.z * gridCellSizes.z) + (gridCellSizes.z/2f)
        );
    }

    // Get the world position of a grid cell based on another (ex. boid's) world position, given bounds an a global grid cell size
    public static Vector3 GetGridCellWorldPositionFromGivenPosition(Vector3Int dimensions, Vector3 origin, float gridCellSize, Vector3 position) {
        Vector3Int xyz = GetGridXYZIndices(dimensions, origin, gridCellSize, position);
        return GetGridCellWorldPositionFromXYZIndices(dimensions, origin, gridCellSize, xyz);
    }

    // Get the world position of a grid cell based on another (ex. boid's) world position, given bounds an a variable grid cell size
    public static Vector3 GetGridCellWorldPositionFromGivenPosition(Vector3Int dimensions, Vector3 origin, Vector3 gridCellSizes, Vector3 position) {
        Vector3Int xyz = GetGridXYZIndices(dimensions, origin, gridCellSizes, position);
        return GetGridCellWorldPositionFromXYZIndices(dimensions, origin, gridCellSizes, xyz);
    }

}

