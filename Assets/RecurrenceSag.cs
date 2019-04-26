using System.Collections.Generic;
using UnityEngine;
using System.Linq;
using System;

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class RecurrenceSag : MonoBehaviour
{
    // User input
    public int subdivisions = 6;
    public float radius = 20.0f;
    public float ClearAperture = 4.0f;
    public float K = 0.0f;
    public Vector3[] MN_coef = new Vector3[1];
    public Vector2[] j_coef = new Vector2[1];
    public Vector3[] MN_calc = new Vector3[1];
    public Boolean otherside;


    public void Awake()
    {
        if (j_coef.Length != 0)
        {
            MN_calc = mnCalc(j_coef);
            GetComponent<MeshFilter>().mesh = RecurrenceSagCreator.Create(subdivisions, radius, ClearAperture, K, MN_calc, otherside);
        }
        else
        {
            GetComponent<MeshFilter>().mesh = RecurrenceSagCreator.Create(subdivisions, radius, ClearAperture, K, MN_coef, otherside);
        }
    }

    private static Vector3[] mnCalc(Vector2[] j_coef)
    {
        Vector3[] mn = new Vector3[j_coef.Length];
        int squareRoot;
        int m_tilda;
        int n_tilda;
        for (int i = 0; i < j_coef.Length; i++)
        {
            squareRoot = (int)Mathf.Sqrt(j_coef[i].x - 1);
            if (((j_coef[i].x - Mathf.Pow(squareRoot, 2)) % 2 != 0))  // odd, cosine Zernike, m>0
            {
                n_tilda = (int)(j_coef[i].x - Mathf.Pow(squareRoot, 2) - 1) / 2;
                m_tilda = squareRoot - n_tilda;
                mn[i].x = m_tilda;
            }
            else
            {
                n_tilda = (int)(j_coef[i].x - Mathf.Pow(squareRoot, 2) - 2) / 2;
                m_tilda = squareRoot - n_tilda;
                mn[i].x = -m_tilda;
            }
            mn[i].y = 2 * n_tilda + m_tilda;
            mn[i].z = j_coef[i].y;
        }
        return mn;
    }

    public void AssignDefaultShader()
    {
        MeshRenderer meshRenderer = (MeshRenderer)gameObject.GetComponent<MeshRenderer>();
        meshRenderer.sharedMaterial = new Material(Shader.Find("Diffuse"));
        meshRenderer.sharedMaterial.color = Color.white;
    }

}

public static class RecurrenceSagCreator
{
    public static Mesh Create(int subdivisions, float radius, float ClearAperture, float K, Vector3[] MN, Boolean otherside)
    {
        //Sanity check for radius and clear aperture
        if (radius < 0) { radius = -radius; }

        //Sanity check for subdivision
        if (subdivisions < 0) { subdivisions = 0; }
        if (subdivisions > 6) { subdivisions = 6; }

        int resolution = 1 << subdivisions;   //2^subdivisions

        //This is set for hemisphere
        Vector3[] vertices = new Vector3[(resolution + 1) * (resolution + 1) * 2 - resolution + 2];
        //Debug.Log ("Total number of vertices as set by resolution: "+vertices.Length.ToString());
        int[] triangles = new int[(1 << (subdivisions * 2 + 2)) * 3];

        // Create Octahedron and divide triangles
        CreateOctahedron(vertices, triangles, resolution, ClearAperture, radius);

        Vector3[] normals = new Vector3[vertices.Length];
        Normalize(vertices, normals, ClearAperture, radius, resolution);

        // Tell if radius is infinity, here we assume when radius > 3000, radius is infinity
        if (radius > 3000)
        {
            for (int i = 0; i < vertices.Length; i++)
            {
                vertices[i].y = 0;
            }
        }

        // Divide m into different categories
        float temp;
        List<List<Vector3>> categoryParent = new List<List<Vector3>>();
        List<Vector3> Input = new List<Vector3>();

        // Convert the input from Vector3 to List
        for (int i = 0; i < MN.Length; i++)
        {
            Input.Add(MN[i]);
        }

        Input.Sort(SortbyX);

        // Put user input of m and n into different categories according to m
        for (int i = 0; i < Input.Count; i++)
        {
            List<Vector3> Category = new List<Vector3>();
            temp = Input[i].x;
            Category.Add(Input[i]);
            for (int j = i + 1; j < Input.Count; j++)
            {
                if (temp == Input[j].x)
                {
                    Category.Add(Input[j]);
                    Input.RemoveAt(j);
                    j--;
                }
            }
            categoryParent.Add(Category);
        }
        //Debug.Log("The first category is:" + categoryParent[0][0]);
        //Debug.Log("The second category is:" + categoryParent[1][0]);

        // Create a new vertex array in order to store the sag departure
        Vector3[] Zernike_vertices = new Vector3[vertices.Length];
        Vector2[] Polar = new Vector2[vertices.Length];

        // System.Array.Copy is important.
        // If we use Array1=Array2, these two arrays will syc all the time because array is only a reference type
        System.Array.Copy(vertices, Zernike_vertices, vertices.Length);

        // Calculate the sag recursively
        for (int j = 0; j < Zernike_vertices.Length; j++)
        {
            Zernike_vertices[j].y = 0;
        }
        Polar = Cart2Pol(Zernike_vertices);
        Zernike_vertices = SagRecurrence(categoryParent, Polar);

        // Add up the zernike y direction departure and calculate the rms over number of points
        float sumZernike = 0;
        for (int i = 0; i < Zernike_vertices.Length; i++)
        {
            sumZernike += Zernike_vertices[i].y;
            //Debug.Log("In recurrence sag, Zernike vertex is: " + Zernike_vertices[i].y);
        }
        //Debug.Log("The sum of y departure is: " + sumZernike);
        sumZernike = Mathf.Sqrt(sumZernike) / Zernike_vertices.Length;
        Debug.Log("In recurrence sag algorithm, the rms of y departure is: " + sumZernike);


        // Scale the surface by radius
        if (radius != 1f)
        {
            for (int i = 0; i < vertices.Length; i++)
            {
                vertices[i] *= radius;
                Zernike_vertices[i] *= radius;
            }
        }

        // Add the base conic
        AddConicTerm(vertices, Zernike_vertices, ClearAperture, radius, K);

        Mesh mesh = new Mesh();
        mesh.name = "Sag Recurrence Zernike";
        mesh.vertices = vertices;
        mesh.normals = normals;
        mesh.triangles = triangles;
        //for (int i = 0; i < normals.Length; i++)
        //    normals[i] = -normals[i];
        mesh.normals = normals;
        if (otherside)
        {
            for (int m = 0; m < mesh.subMeshCount; m++)
            {
                for (int i = 0; i < triangles.Length; i += 3)
                {
                    int t = triangles[i + 0];
                    triangles[i + 0] = triangles[i + 1];
                    triangles[i + 1] = t;
                }
                mesh.SetTriangles(triangles, m);
            }
        }
        return mesh;
    }

    private static int SortbyX(Vector3 x, Vector3 y)
    {
        return x.x.CompareTo(y.x);
    }

    private static void CreateOctahedron(Vector3[] vertices, int[] triangles, int resolution, float ClearAperture, float Radius)
    {
        int v = 0, vBottom = 0, t = 0;           //v tracks the number of vertices
                                                 //vBottom tracks the first vertex index of the previous row

        for (int i = 0; i < 4; i++)
        {
            vertices[v++] = Vector3.down;
        }
        for (int i = 1; i <= resolution; i++)
        {
            float progress = (float)i / (resolution);
            //Debug.Log ("Show progress: "+progress);
            Vector3 from, to;
            //Stupid method
            Vector3 front = new Vector3(0, ClearAperture / (2 * Radius) - 1, ClearAperture / (2 * Radius));
            Vector3 left = new Vector3(-ClearAperture / (2 * Radius), ClearAperture / (2 * Radius) - 1, 0);
            Vector3 back = new Vector3(0, ClearAperture / (2 * Radius) - 1, -ClearAperture / (2 * Radius));
            Vector3 right = new Vector3(ClearAperture / (2 * Radius), ClearAperture / (2 * Radius) - 1, 0);
            Vector3 forward = new Vector3(0, ClearAperture / (2 * Radius) - 1, ClearAperture / (2 * Radius));

            vertices[v++] = to = Vector3.Lerp(Vector3.down, front, progress);

            from = to;
            to = Vector3.Lerp(Vector3.down, left, progress);
            t = CreateLowerStrip(i, v, vBottom, t, triangles);
            v = CreateVertexLine(from, to, i, v, vertices);
            vBottom += i > 1 ? (i - 1) : 1;

            from = to;
            to = Vector3.Lerp(Vector3.down, back, progress);
            t = CreateLowerStrip(i, v, vBottom, t, triangles);
            v = CreateVertexLine(from, to, i, v, vertices);
            vBottom += i > 1 ? (i - 1) : 1;

            from = to;
            to = Vector3.Lerp(Vector3.down, right, progress);
            t = CreateLowerStrip(i, v, vBottom, t, triangles);
            v = CreateVertexLine(from, to, i, v, vertices);
            vBottom += i > 1 ? (i - 1) : 1;

            from = to;
            to = Vector3.Lerp(Vector3.down, forward, progress);
            t = CreateLowerStrip(i, v, vBottom, t, triangles);
            v = CreateVertexLine(from, to, i, v, vertices);
            vBottom += i > 1 ? (i - 1) : 1;

            vBottom = v - 1 - i * 4;
        }
    }

    private static int CreateLowerStrip(int steps, int vTop, int vBottom, int t, int[] triangles)
    {
        for (int i = 1; i < steps; i++)
        {
            triangles[t++] = vBottom;
            triangles[t++] = vTop;
            triangles[t++] = vTop - 1;

            triangles[t++] = vBottom++;
            triangles[t++] = vBottom;
            triangles[t++] = vTop++;
        }

        triangles[t++] = vBottom;
        triangles[t++] = vTop;
        triangles[t++] = vTop - 1;
        return t;
    }


    private static int CreateVertexLine(
        Vector3 from, Vector3 to, int steps, int v, Vector3[] vertices
    )
    {
        for (int i = 1; i <= steps; i++)
        {
            vertices[v++] = Vector3.Lerp(from, to, (float)i / steps);
        }
        return v;
    }


    private static void Normalize(Vector3[] vertices, Vector3[] normals, float ClearAperture, float radius, int resolution)
    {
        float Ratio; // The ratio of each vertex y value over the last ring y value
        float totalHeight; // Last ring y value + 1 (MP)
        float tempHeight; // totalHeight * Ratio
        float tempYvalue; // tempHeight + last ring y value
        float RingRadius;
        float One;
        Vector3 lastringYvalueBefore;
        Vector3 lastringYvalueAfter;

        // Before normalization, find the last ring vertices y value
        lastringYvalueBefore = vertices[0];
        for (int i = 0; i < vertices.Length; i++)
        {
            if (vertices[i].y > lastringYvalueBefore.y)
            {

                lastringYvalueBefore = vertices[i];
            }
        }
        //Debug.Log ("Before normalization, last ring y value is: " + lastringYvalueBefore.y);

        // Find the last ring y value after normalization
        Vector3 temp1 = vertices[vertices.Length - 1];
        // Find the last ring y value after normalization
        RingRadius = ClearAperture / (2 * radius);
        One = Mathf.Sqrt(temp1.x * temp1.x + temp1.z * temp1.z);
        temp1.x = temp1.x / One * RingRadius;
        temp1.z = temp1.z / One * RingRadius;
        temp1.y = -Mathf.Sqrt(1 - temp1.x * temp1.x - temp1.z * temp1.z);
        lastringYvalueAfter = temp1;
        //Debug.Log ("After normalization, last ring y value is: " + lastringYvalueAfter.y);

        for (int i = 0; i < vertices.Length; i++)
        {
            Vector3 temp = vertices[i];
            //Debug.Log ("Before normalization, vertex is: " + temp);
            // If vertex coordinate y value is -1, just normalize it
            if (temp.y == -1)
            {
                normals[i] = vertices[i] = vertices[i].normalized;
            }
            // If vertex coordinate y value equals to last ring vertex y value, normalize it to the ring
            else if (temp.y == lastringYvalueBefore.y)
            {
                // Find the last ring y value after normalization
                RingRadius = ClearAperture / (2 * radius);
                One = Mathf.Sqrt(temp.x * temp.x + temp.z * temp.z);
                temp.x = temp.x / One * RingRadius;
                temp.z = temp.z / One * RingRadius;
                temp.y = -Mathf.Sqrt(1 - temp.x * temp.x - temp.z * temp.z);
                if (float.IsNaN(temp.y))
                {
                    temp.y = 0;
                }
                normals[i] = vertices[i] = temp;
                //Debug.Log ("The vertices are: " +temp);
            }
            else
            {
                Ratio = (temp.y + 1) / (lastringYvalueBefore.y + 1);
                totalHeight = lastringYvalueAfter.y + 1;
                tempHeight = totalHeight * Ratio;
                tempYvalue = lastringYvalueAfter.y + tempHeight - totalHeight; // last ring y value is a negative value
                RingRadius = Mathf.Sqrt(1 - Mathf.Pow(tempYvalue, 2));
                One = Mathf.Sqrt(temp.x * temp.x + temp.z * temp.z);
                temp.x = temp.x / One * RingRadius;
                temp.z = temp.z / One * RingRadius;
                temp.y = tempYvalue;
                normals[i] = vertices[i] = temp;
            }

        }
    }

    private static Vector2[] Cart2Pol(Vector3[] vertices)
    {
        Vector2[] Polar = new Vector2[vertices.Length];
        for (int i = 0; i < vertices.Length; i++)
        {
            Polar[i].x = Mathf.Sqrt(Mathf.Pow(vertices[i].x, 2) + Mathf.Pow(vertices[i].z, 2));
            Polar[i].y = Mathf.Atan2(vertices[i].z, vertices[i].x);
            // Debug.Log ("The Polar coordinates are: " +Polar [i].x);
        }
        return Polar;
    }


    private static float[,] MatrixMultiplication(float[,] Matrix1, float[,] Matrix2)
    {
        float[,] Result = new float[Matrix1.GetLength(0), Matrix2.GetLength(1)];
        if (Matrix1.GetLength(1) == Matrix2.GetLength(0))
        {
            for (int i = 0; i < Result.GetLength(0); i++)
            {
                for (int j = 0; j < Result.GetLength(1); j++)
                {
                    Result[i, j] = 0;
                    for (int k = 0; k < Matrix1.GetLength(1); k++)
                    {
                        Result[i, j] = Result[i, j] + Matrix1[i, k] * Matrix2[k, j];
                    }
                }
            }
            return Result;
        }
        else
        {
            for (int i = 0; i < Result.GetLength(0); i++)
            {
                for (int j = 0; j < Result.GetLength(1); j++)
                {
                    Result[i, j] = 0;
                }
            }
            return Result;
        }
    }

    private static Vector3[] SagRecurrence(List<List<Vector3>> Category, Vector2[] Polar)
    { // Each Category[i][j] is a Vector3 (m,n,c), each Category [i] has same m value, sorted from smaller m to larger m
        int M = Category.Count;   // Number of category, starts at 0
        float[] Rho = new float[Polar.Length];
        float[] theta = new float[Polar.Length];

        float[] zero = new float[Polar.Length];
        List<float[]> m_vertices = new List<float[]>();

        Vector3[] S = new Vector3[Polar.Length];
        Vector3[] temp_cos = new Vector3[Polar.Length];
        Vector3[] temp_sin = new Vector3[Polar.Length];

        for (int i = 0; i < Polar.Length; i++)
        {
            Rho[i] = Polar[i].x;
            theta[i] = Polar[i].y;
            zero[i] = 0;
        }


        // Check if the first input for m == 0, if m==0, Do it outside the sum
        // Category is sorted so that only the category[0] term can have m=0
        int added = 0;
        if (Category[0][0].x == 0)
        {
            Debug.Log("first term is 0");
            m_vertices.Add(Clenshaw(Category[0], Rho));
            //Debug.Log("The polar length is:" + Polar.Length);
            //for (int j = 0; j < Polar.Length; j++)
            //{

            //    Debug.Log("The m_vertices at the beginning is: " + m_vertices[0][j]);
            //}
        }
        else if (Category[0][0].x != 0)
        {
            m_vertices.Add(zero); added++;
            //Debug.Log("The polar length is:" + Polar.Length);
            m_vertices.Add(Clenshaw(Category[0], Rho));
            Debug.Log("first term is not 0, added zero");
        }
        for (int i = 1; i < M; i++)
        {
            m_vertices.Add(Clenshaw(Category[i], Rho));
            Debug.Log("added term");
            //for (int j = 0; j < Polar.Length; j++)
            //{
            //Debug.Log("When recurrence sag, the m_vertices at the beginning is: " + m_vertices[1][j]);
            //}
        }
        // calcualte the sag
        for (int i = 0; i < Polar.Length; i++)
        {
            for (int j = 0; j < M; j++)
            {
                if (Category[j][0].x >= 0)
                {
                    temp_cos[i].y += m_vertices[j + added][i] * Mathf.Pow(Rho[i], Category[j][0].x) * Mathf.Cos(theta[i] * Mathf.Abs(Category[j][0].x));
                    //Debug.Log("m_vertices is: " + m_vertices[j][i]);
                    //Debug.Log("temp_cos is: " + temp_cos[i].y);
                }
                else if (Category[j][0].x < 0)
                {
                    temp_sin[i].y += m_vertices[j + added][i] * Mathf.Pow(Rho[i], Mathf.Abs(Category[j][0].x)) * Mathf.Sin(theta[i] * Mathf.Abs(Category[j][0].x));
                }
            }
            S[i].y = m_vertices[0][i] + temp_cos[i].y + temp_sin[i].y;
        }

        return S;
    }

    private static float[] Clenshaw(List<Vector3> MN, float[] rho)
    {
        float m_tilda = Mathf.Abs(MN[0].x); //MN is one row in category, having the same m value
        int max_y = int.MinValue;         // find the max of n to calculate N
        float max_z = float.MinValue;
        for (int i = 0; i < MN.Count; i++)
        {
            if ((int)MN[i].y > max_y)
            {
                max_y = (int)MN[i].y; // max n
                max_z = MN[i].z;
            }
        }

        int N = (int)(max_y - m_tilda) / 2; //Maximum n_tilda
        int[] n_tilda = new int[MN.Count];
        float[] coefficient = new float[N + 1];
        //Debug.Log("The mn count is: " + MN.Count);
        //Debug.Log("N is: " + N);

        int N_1 = N - 1;  // Record the n for N-1 term

        //calculate n_tilda
        for (int i = 0; i < MN.Count; i++)
        {
            n_tilda[i] = (int)(MN[i].y - m_tilda) / 2;
        }

        // get an array for coefficient
        for (int i = N; i >= 0; i--)
        {
            for (int j = 0; j < n_tilda.Length; j++)
            {
                if (n_tilda[j] == i)
                {
                    coefficient[i] += MN[j].z;
                }
                else
                {
                    coefficient[i] += 0;
                }
            }
            Debug.Log("The coefficient is: " + coefficient[i]);
        }

        float a_N1;
        float b_N1;
        float c_N1;
        float s_N1;      // Record the s for N-1 term

        float a;
        float b;
        float c;  // actually c calculate c+1 in the loop
        float s;
        float s_C1; // to calculate c_n+1

        // calculate the a1 b1 c1 c2
        float a1;
        float b1;
        float c1;
        float c2;
        float s1;
        float s2;

        float[] Sag = new float[rho.Length];
        float[,] alpha = new float[rho.Length, N + 1];

        if (N == 0)
        {
            for (int i = 0; i < rho.Length; i++)
            {
                // Add up all the (0,0) terms coefficient
                for (int j = 0; j < MN.Count; j++)
                {
                    Sag[i] += MN[j].z;
                }
            }
        }

        if (N == 1)
        {
            // there are two terms here, the lower order in y is the first term and the higher order in y is the second term
            float y_lower = 0;
            float y_higher = 0;
            int min = int.MinValue;
            int max = int.MinValue;
            for (int j = 0; j < MN.Count; j++)
            {
                if (min > MN[j].y) { min = (int)MN[j].y; }
                if (max < MN[j].y) { max = (int)MN[j].y; }

                if (MN[j].y == min)
                {
                    y_lower += MN[j].z;
                }
                if (MN[j].y == max)
                {
                    y_higher += MN[j].z;
                }
            }

            for (int i = 0; i < rho.Length; i++)
            {
                Sag[i] = y_lower + y_higher * (2 * Mathf.Pow(rho[i], 2) - m_tilda + m_tilda * Mathf.Pow(rho[i], 2) - 1);
                //Debug.Log("Sag is: " + Sag[i]);
            }
        }

        if (N == 2)
        {
            float y_lower = 0;
            float y_middle = 0;
            float y_higher = 0;
            int min = int.MinValue;
            int max = int.MinValue;
            for (int j = 0; j < MN.Count; j++)
            {
                if (min > MN[j].y) { min = (int)MN[j].y; }
                if (max < MN[j].y) { max = (int)MN[j].y; }
                if (MN[j].y == min)
                {
                    y_lower += MN[j].z;
                    //Debug.Log("y_lower is: " + y_lower);
                }

                else if (MN[j].y == max)
                {
                    y_higher += MN[j].z;
                    //Debug.Log("y_higher is: " + y_higher);
                }
                else
                {
                    y_middle += MN[j].z;
                    //Debug.Log("y_middle is: " + y_middle);
                }
            }
            for (int i = 0; i < rho.Length; i++)
            {
                Sag[i] = y_lower + y_middle * (2 * Mathf.Pow(rho[i], 2) - m_tilda + m_tilda * Mathf.Pow(rho[i], 2) - 1) + y_higher * ((Mathf.Pow(m_tilda, 2) / 2 + 7 * m_tilda / 2 + 6) * Mathf.Pow(rho[i], 4) - (Mathf.Pow(m_tilda, 2) + 5 * m_tilda + 6) * Mathf.Pow(rho[i], 2) + Mathf.Pow(m_tilda, 2) / 2 + 3 * m_tilda / 2 + 1);
                //Debug.Log("In recurrence sag, the sag is:" + Sag[i]);
            }
        }

        if (N > 2)
        {
            s_N1 = m_tilda + 2 * N_1;
            a_N1 = -(s_N1 + 1) * (Mathf.Pow((s_N1 - N_1), 2) + Mathf.Pow(N_1, 2) + s_N1) / ((N_1 + 1) * (s_N1 - N_1 + 1) * s_N1);
            b_N1 = (s_N1 + 2) * (s_N1 + 1) / ((N_1 + 1) * (s_N1 - N_1 + 1));
            c_N1 = (s_N1 + 2) * (s_N1 - N_1) * N_1 / ((N_1 + 1) * (s_N1 - N_1 + 1) * s_N1);
            Debug.Log("a_N1 is: " + a_N1);
            Debug.Log("b_N1 is: " + b_N1);
            Debug.Log("c_N1 is: " + c_N1);

            s1 = m_tilda + 2;
            s2 = m_tilda + 4;
            a1 = -(s1 + 1) * (Mathf.Pow((s1 - 1), 2) + 1 + s1) / (2 * s1 * s1);
            b1 = (s1 + 2) * (s1 + 1) / (2 * s1);
            c1 = (s1 + 2) * (s1 - 1) / (2 * s1 * s1);
            c2 = (s2 + 2) * (s2 - 2) * 2 / (3 * (s2 - 1) * s2);
            Debug.Log("a1 is: " + a1);
            Debug.Log("b1 is: " + b1);
            Debug.Log("c1 is: " + c1);
            Debug.Log("c2 is: " + c2);

            for (int i = 0; i < rho.Length; i++)
            {
                // Define the first two terms of alpha
                alpha[i, N] = max_z;
                alpha[i, N - 1] = coefficient[N - 1] + (a_N1 + b_N1 * Mathf.Pow(rho[i], 2)) * max_z;
                //Debug.Log("alpha(N) is: " + alpha[i, N]);
                //Debug.Log("alpha(N-1) is: " + alpha[i, N-1]);
            }
            if (N > 3)
            {
                for (int j = N - 2; j >= 2; j--)
                {
                    s = m_tilda + 2 * j;
                    s_C1 = m_tilda + 2 * (j + 1);
                    a = -(s + 1) * (Mathf.Pow((s - j), 2) + Mathf.Pow(j, 2) + s) / ((j + 1) * (s - j + 1) * s);
                    b = (s + 2) * (s + 1) / ((j + 1) * (s - j + 1));
                    c = (s_C1 + 2) * (s_C1 - (j + 1)) * (j + 1) / (((j + 1) + 1) * (s_C1 - (j + 1) + 1) * s_C1);
                    //Debug.Log("a is: " + a);
                    //Debug.Log("b is: " + b);
                    //Debug.Log("c is: " + c);
                    //Debug.Log("The coefficient for N>3 is: " + coefficient[j]);
                    for (int i = 0; i < rho.Length; i++)
                    {
                        alpha[i, j] = coefficient[j] + (a + b * Mathf.Pow(rho[i], 2)) * alpha[i, j + 1] - c * alpha[i, j + 2];
                        //Debug.Log("alpha4 is: " + alpha[i, j+2]);
                        //Debug.Log("alpha3 is: " + alpha[i, j+1]);
                        //Debug.Log("alpha2 is: " + alpha[i, j]);
                    }
                }
            }
            float[] aa = new float[rho.Length];
            for (int i = 0; i < rho.Length; i++)
            {
                Sag[i] = (alpha[i, 2] * (a1 + b1 * Mathf.Pow(rho[i], 2)) + coefficient[1] - c2 * alpha[i, 3]) * (2 * Mathf.Pow(rho[i], 2) - m_tilda + m_tilda * Mathf.Pow(rho[i], 2) - 1) + (coefficient[0] - alpha[i, 2] * c1);
                //Debug.Log("Sag is: " + Sag[i]);
            }
        }
        return Sag;
    }

    private static void AddConicTerm(Vector3[] vertices, Vector3[] Zernike_vertices, float ClearAperture, float radius, float K)
    {
        float z; // z is the variable in the equation

        for (int i = 0; i < vertices.Length; i++)
        {
            Vector3 temp = vertices[i];
            float r = Mathf.Sqrt(temp.x * temp.x + temp.z * temp.z);
            z = Mathf.Pow(r, 2) / (radius * (1 + Mathf.Sqrt(1 - (1 + K) * Mathf.Pow(r, 2) / (Mathf.Pow(radius, 2)))));
            // Check if z is NAN because of the denominator
            if (System.Single.IsNaN(z))
            {
                z = radius;
            }


            //Debug.Log ("The Zernike term y coordinate is" + Zernike_vertices[i].y);
            temp.y = -z + radius + Zernike_vertices[i].y;
            vertices[i] = temp;
            //Debug.Log ("The sum y coordinate is" + temp.y);
        }
    }
}
