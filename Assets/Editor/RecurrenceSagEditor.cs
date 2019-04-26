using UnityEngine;
using UnityEditor;
using System.Collections;

[CustomEditor(typeof(RecurrenceSag))]

public class RecurrenceSagEditor : Editor
{
    [MenuItem("GameObject/Create Other/Recurrence Sag")]
    static void Create()
    {
        GameObject gameObject = new GameObject("RecurrenceSag");
        RecurrenceSag c = gameObject.AddComponent<RecurrenceSag>();
        MeshFilter meshFilter = gameObject.GetComponent<MeshFilter>();
        meshFilter.mesh = new Mesh();
        c.AssignDefaultShader();
        c.Awake();
    }

    public override void OnInspectorGUI()
    {
        RecurrenceSag obj;
        obj = target as RecurrenceSag;
        if (obj == null)
        {
            return;
        }

        base.DrawDefaultInspector();
        //if (GUILayout.Button("Standard Zernike"))
        //{

        //}
        //if (GUILayout.Button("Fringe Zernike"))
        //{
        //    //add everthing the button would do.
        //}
        if (GUI.changed)
        {
            obj.Awake();
        }
    }
}