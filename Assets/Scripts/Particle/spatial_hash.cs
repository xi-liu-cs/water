using UnityEngine;

public class spatial_hash
{
    public static Vector3Int dimension;
    public static float grid_size;
    public static float[] bound;

    public static Vector3Int get_cell(Vector3 position)
    {
        return new Vector3Int((int)((position.x - bound[0]) / grid_size), (int)((position.y - bound[2]) / grid_size), (int)((position.z - bound[4]) / grid_size));
    }

    public static int hash(Vector3Int cell)
    {
        return cell.x + dimension.x * (cell.y + dimension.y * cell.z);
    }
}