using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class StoredValues : MonoBehaviour {
    public float[] zernike;
    public List<float> oldValues;
    public RecurrenceSag side1, side2;
    public Material green, blue; // (Prefab) The material of the selected obj
    public Material white; // (Prefab) Used during selection to turn the unselected obj back to white
    public Vector3 pos;
	// Use this for initialization
	void Start () {
        zernike = new float[37];
        oldValues = new List<float>(37);
    }
	
	// Update is called once per frame
	void Update () {
		
	}
    public void edit(int term, float value)
    {
        zernike[term] = value;
        load();
    }
    public void confirm()
    {
        print("confirm");
        oldValues = new List<float>();
        for(int i = 0; i < 37; i++) { oldValues.Add(0); }
        for (int i = 0; i < zernike.Length; i++) oldValues[i]=zernike[i]; 
    }
    public void undo()
    {
        print("undo");
        zernike = new float[37];
        for (int i = 0; i < oldValues.Count; i++)
            zernike[i] = oldValues[i];
        load();
    }
    public void empty()
    {
        zernike = new float[37];
        load();
    }
    public void load() {
        side1.transform.localPosition = new Vector3(0, -side1.radius, 0);
        side2.transform.localPosition = new Vector3(0, -side1.radius, 0);
        int n = 0; //the number of non-zero terms
        // j[] stores the non-zero terms
        // Since length of a list[] can not be changed, need to determine its length first
        for(int i = 0; i < zernike.Length; i++)
        {
            if(zernike[i] != 0) n++;
        }
        Vector2[] j = new Vector2[n]; // Create empty j[] with n terms
        n = 0;
        // Now add numbers to j
        for (int i = 0; i < zernike.Length; i++)
        {
            if (zernike[i] != 0) { n++;
                if (i == 4) {
                    side1.transform.localPosition += new Vector3(0, side1.radius * 2 * zernike[i], 0);
                    side2.transform.localPosition += new Vector3(0, side1.radius * 2 * zernike[i], 0);
                }
                j[n - 1] = new Vector2(i, zernike[i]); 
            }
        }
        side1.j_coef = j;
        side1.Awake();
        side2.j_coef = j;
        side2.Awake();
        
    }
    public float[] getZernike() {
        return zernike;
    }
    public void setWhite() {
        side1.GetComponent<Renderer>().material = white;
        side2.GetComponent<Renderer>().material = white;
    }
    public void setGreen() {
        side1.GetComponent<Renderer>().material = green;
        side2.GetComponent<Renderer>().material = green;
    }
    public void setBlue() {
        side1.GetComponent<Renderer>().material = blue;
        side2.GetComponent<Renderer>().material = blue;
    }
}
