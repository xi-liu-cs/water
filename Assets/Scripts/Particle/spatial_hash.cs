using UnityEngine;

public class spatial_hash
{
    public static Vector3Int dimension;
    public static float grid_size;

    public static Vector3Int get_cell(Vector3 position)
    {
        return new Vector3Int((int)(position.x / grid_size), (int)(position.y / grid_size), (int)(position.z / grid_size));
    }

    public static int hash(Vector3Int cell)
    {
        return cell.x + dimension.x * (cell.y + dimension.y * cell.z);
    }
}