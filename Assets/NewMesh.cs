using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class NewMesh : MonoBehaviour
{
   // Create a new Mesh and assign it to the Mesh filter
    void Start()
    {
        Mesh mesh = new Mesh();
        GetComponent<MeshFilter>().mesh = mesh;
    }
}
