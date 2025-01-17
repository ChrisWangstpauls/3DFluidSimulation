using UnityEngine;
using System.Collections.Generic;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;
using Unity.Burst;

public class FluidSimulation : MonoBehaviour
{
	[Header("Simulation Parameters")]
	public int size = 100;
	public float diffusion = 0.0001f;
	public float viscosity = 0.0001f;
	public float timeStep = 0.1f;

	[Header("Visualization")]
	public Color fluidColor = Color.white;
	[Range(0f, 1f)]
	public float colorIntensity = 1f;
	public bool useGradient = false;
	public Gradient colorGradient;

	private float[] density;
	private float[] velocityX;
	private float[] velocityY;
	private float[] velocityX0;
	private float[] velocityY0;

	private Material fluidMaterial;
	private Texture2D fluidTexture;
	private Camera mainCamera;
	private Vector3[] quadCorners = new Vector3[4];

	void OnValidate()
	{
		// Update visualization when parameters change in editor
		if (fluidTexture != null)
		{
			UpdateVisualization();
		}
	}

	void Start()
	{
		// Initialize arrays
		int totalSize = size * size;
		density = new float[totalSize];
		velocityX = new float[totalSize];
		velocityY = new float[totalSize];
		velocityX0 = new float[totalSize];
		velocityY0 = new float[totalSize];

		// Create visualization texture
		fluidTexture = new Texture2D(size, size);
		fluidTexture.filterMode = FilterMode.Point;

		// Create material for rendering
		fluidMaterial = new Material(Shader.Find("Unlit/Texture"));
		fluidMaterial.mainTexture = fluidTexture;

		// Create quad for visualization
		GameObject quad = GameObject.CreatePrimitive(PrimitiveType.Quad);
		quad.GetComponent<Renderer>().material = fluidMaterial;

		mainCamera = Camera.main;

		// Store quad mesh filter for corner calculations
		MeshFilter meshFilter = quad.GetComponent<MeshFilter>();
		List<Vector3> cornersList = new List<Vector3>();
		meshFilter.mesh.GetVertices(cornersList);
		quadCorners = cornersList.ToArray();

		// Transform corners to world space
		for (int i = 0; i < 4; i++)
		{
			quadCorners[i] = quad.transform.TransformPoint(quadCorners[i]);
		}

		// Initialize default gradient if none is set
		if (colorGradient.Equals(new Gradient()))
		{
			GradientColorKey[] colorKeys = new GradientColorKey[2];
			colorKeys[0].color = Color.blue;
			colorKeys[0].time = 0.0f;
			colorKeys[1].color = Color.red;
			colorKeys[1].time = 1.0f;

			GradientAlphaKey[] alphaKeys = new GradientAlphaKey[2];
			alphaKeys[0].alpha = 1.0f;
			alphaKeys[0].time = 0.0f;
			alphaKeys[1].alpha = 1.0f;
			alphaKeys[1].time = 1.0f;

			colorGradient.SetKeys(colorKeys, alphaKeys);
		}
	}

	void Update()
	{
		// Add forces based on mouse input
		if (Input.GetMouseButton(0))
		{
			Vector2 mousePos = GetMousePositionInGrid();
			if (mousePos.x >= 0 && mousePos.x < size && mousePos.y >= 0 && mousePos.y < size)
			{
				AddDensity(mousePos.x, mousePos.y, 100);
				AddVelocity(mousePos.x, mousePos.y, Input.GetAxis("Mouse X") * 10, Input.GetAxis("Mouse Y") * 10);
			}
		}

		Simulate();
		UpdateVisualization();
	}

	Vector2 GetMousePositionInGrid()
	{
		Vector3 mouseScreenPos = Input.mousePosition;
		Vector3 mouseWorldPos = mainCamera.ScreenToWorldPoint(new Vector3(mouseScreenPos.x, mouseScreenPos.y, -mainCamera.transform.position.z));

		float minX = Mathf.Min(quadCorners[0].x, quadCorners[1].x, quadCorners[2].x, quadCorners[3].x);
		float maxX = Mathf.Max(quadCorners[0].x, quadCorners[1].x, quadCorners[2].x, quadCorners[3].x);
		float minY = Mathf.Min(quadCorners[0].y, quadCorners[1].y, quadCorners[2].y, quadCorners[3].y);
		float maxY = Mathf.Max(quadCorners[0].y, quadCorners[1].y, quadCorners[2].y, quadCorners[3].y);

		float normalizedX = (mouseWorldPos.x - minX) / (maxX - minX);
		float normalizedY = (mouseWorldPos.y - minY) / (maxY - minY);

		return new Vector2(normalizedX * size, normalizedY * size);
	}

	void Simulate()
	{
		VelocityStep();
		DensityStep();
	}

	void VelocityStep()
	{
		Diffuse(1, velocityX0, velocityX, viscosity);
		Diffuse(2, velocityY0, velocityY, viscosity);

		Project(velocityX0, velocityY0, velocityX, velocityY);

		Advect(1, velocityX, velocityX0, velocityX0, velocityY0);
		Advect(2, velocityY, velocityY0, velocityX0, velocityY0);

		Project(velocityX, velocityY, velocityX0, velocityY0);
	}

	void DensityStep()
	{
		float[] densityTemp = new float[size * size];
		Diffuse(0, densityTemp, density, diffusion);
		Advect(0, density, densityTemp, velocityX, velocityY);
	}

	void AddDensity(float x, float y, float amount)
	{
		int i = Mathf.Clamp((int)x, 0, size - 1);
		int j = Mathf.Clamp((int)y, 0, size - 1);

		density[IX(i, j)] += amount;
	}

	void AddVelocity(float x, float y, float amountX, float amountY)
	{
		int i = Mathf.Clamp((int)x, 0, size - 1);
		int j = Mathf.Clamp((int)y, 0, size - 1);

		velocityX[IX(i, j)] += amountX;
		velocityY[IX(i, j)] += amountY;
	}

	void Diffuse(int b, float[] x, float[] x0, float diff)
	{
		float a = timeStep * diff * (size - 2) * (size - 2);
		LinearSolve(b, x, x0, a, 1 + 6 * a);
	}

	void LinearSolve(int b, float[] x, float[] x0, float a, float c)
	{
		for (int k = 0; k < 20; k++)
		{
			for (int i = 1; i < size - 1; i++)
			{
				for (int j = 1; j < size - 1; j++)
				{
					x[IX(i, j)] = (x0[IX(i, j)] + a * (
						x[IX(i + 1, j)] + x[IX(i - 1, j)] +
						x[IX(i, j + 1)] + x[IX(i, j - 1)]
					)) / c;
				}
			}
			SetBoundary(b, x);
		}
	}

	void Project(float[] velocX, float[] velocY, float[] p, float[] div)
	{
		for (int i = 1; i < size - 1; i++)
		{
			for (int j = 1; j < size - 1; j++)
			{
				div[IX(i, j)] = -0.5f * (
					velocX[IX(i + 1, j)] - velocX[IX(i - 1, j)] +
					velocY[IX(i, j + 1)] - velocY[IX(i, j - 1)]
				) / size;
				p[IX(i, j)] = 0;
			}
		}

		SetBoundary(0, div);
		SetBoundary(0, p);
		LinearSolve(0, p, div, 1, 6);

		for (int i = 1; i < size - 1; i++)
		{
			for (int j = 1; j < size - 1; j++)
			{
				velocX[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * size;
				velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * size;
			}
		}

		SetBoundary(1, velocX);
		SetBoundary(2, velocY);
	}

	void Advect(int b, float[] d, float[] d0, float[] velocX, float[] velocY)
	{
		float dt0 = timeStep * (size - 2);

		for (int i = 1; i < size - 1; i++)
		{
			for (int j = 1; j < size - 1; j++)
			{
				float x = i - dt0 * velocX[IX(i, j)];
				float y = j - dt0 * velocY[IX(i, j)];

				if (x < 0.5f) x = 0.5f;
				if (x > size - 1.5f) x = size - 1.5f;
				int i0 = (int)x;
				int i1 = i0 + 1;

				if (y < 0.5f) y = 0.5f;
				if (y > size - 1.5f) y = size - 1.5f;
				int j0 = (int)y;
				int j1 = j0 + 1;

				float s1 = x - i0;
				float s0 = 1 - s1;
				float t1 = y - j0;
				float t0 = 1 - t1;

				d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
							  s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
			}
		}

		SetBoundary(b, d);
	}

	void SetBoundary(int b, float[] x)
	{
		for (int i = 1; i < size - 1; i++)
		{
			x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
			x[IX(size - 1, i)] = b == 1 ? -x[IX(size - 2, i)] : x[IX(size - 2, i)];
			x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
			x[IX(i, size - 1)] = b == 2 ? -x[IX(i, size - 2)] : x[IX(i, size - 2)];
		}

		x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
		x[IX(0, size - 1)] = 0.5f * (x[IX(1, size - 1)] + x[IX(0, size - 2)]);
		x[IX(size - 1, 0)] = 0.5f * (x[IX(size - 2, 0)] + x[IX(size - 1, 1)]);
		x[IX(size - 1, size - 1)] = 0.5f * (x[IX(size - 2, size - 1)] + x[IX(size - 1, size - 2)]);
	}

	int IX(int x, int y)
	{
		return x + y * size;
	}

	void UpdateVisualization()
	{
		Color[] colors = new Color[size * size];

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				float d = density[IX(i, j)] * colorIntensity;
				if (useGradient)
				{
					// Use gradient evaluation (d is clamped to 0-1 range)
					colors[IX(i, j)] = colorGradient.Evaluate(Mathf.Clamp01(d));
				}
				else
				{
					// Use single color with varying intensity
					colors[IX(i, j)] = new Color(
						fluidColor.r * d,
						fluidColor.g * d,
						fluidColor.b * d,
						fluidColor.a
					);
				}
			}
		}

		fluidTexture.SetPixels(colors);
		fluidTexture.Apply();
	}
}
