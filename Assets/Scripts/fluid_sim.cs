using UnityEngine;
using System.Collections.Generic;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;
using Unity.Burst;

public class FluidSimulation : MonoBehaviour
{
	[Header("Simulation Parameters")]
	[Range(32, 512)]
	public int size = 128;
	[Tooltip("Physical size of the simulation area")]
	public float physicalSize = 1.0f;
	[Range(0.1f, 10f)]
	public float resolutionMultiplier = 1.0f;
	public float diffusion = 0.0001f;
	public float viscosity = 0.0001f;
	public float timeStep = 0.1f;
	public bool autoAdjustParameters = true;

	[Header("Customizable Source")]
	public bool enableCustomSource = false;
	[Range(1f, 500f)]
	public float sourceStrength = 100f;
	public bool sourceEmitsVelocity = false;
	[Range(0f, 360f)]
	public float sourceDirection = 0f;
	[Range(1f, 50f)]
	public float sourceVelocity = 10f;
	[Range(0.1f, 5f)]
	public float sourceRadius = 1f;
	[Range(0.1f, 5f)]
	public float sourcePulseRate = 1f;
	public bool sourcePulsing = false;
	[Range(0f, 1f)]
	public float sourcePositionX = 0.5f; // Normalized position (0-1)
	[Range(0f, 1f)]
	public float sourcePositionY = 0.5f; // Normalized position (0-1)
	public bool moveSourceWithMouse = false;
	public KeyCode sourcePositionKey = KeyCode.LeftShift; // Hold this key to position source with mouse

	[Header("Visualization")]
	public Color fluidcolour = Color.white;
	[Range(0f, 1f)]
	public float colourIntensity = 1f;
	public bool useGradient = false;
	public Gradient colourGradient;
	public bool visualizeSourcePosition = true;
	public Color sourcePositionColor = Color.yellow;

	private float[] density;
	private float[] velocityX;
	private float[] velocityY;
	private float[] velocityX0;
	private float[] velocityY0;

	private int currentSize;
	private float cellSize;
	private float dtScale;
	private float elapsedTime = 0f;

	private Material fluidMaterial;
	private Texture2D fluidTexture;
	private Camera mainCamera;
	private Vector3[] quadCorners = new Vector3[4];

	private int previousSize;
	private float previousResolutionMultiplier;

	void OnValidate()
	{
		// Update resolution when parameters change in editor
		if (previousSize != size || previousResolutionMultiplier != resolutionMultiplier)
		{
			// Only reinitialize if we've already started
			if (fluidTexture != null)
			{
				ResetSimulation();
			}

			previousSize = size;
			previousResolutionMultiplier = resolutionMultiplier;
		}

		// Update visualization when parameters change in editor
		if (fluidTexture != null)
		{
			UpdateVisualization();
		}
	}

	void Start()
	{
		ResetSimulation();

		// Initialize default gradient if none is set
		if (colourGradient.Equals(new Gradient()))
		{
			GradientColorKey[] colourKeys = new GradientColorKey[2];
			colourKeys[0].color = Color.blue;
			colourKeys[0].time = 0.0f;
			colourKeys[1].color = Color.red;
			colourKeys[1].time = 1.0f;

			GradientAlphaKey[] alphaKeys = new GradientAlphaKey[2];
			alphaKeys[0].alpha = 1.0f;
			alphaKeys[0].time = 0.0f;
			alphaKeys[1].alpha = 1.0f;
			alphaKeys[1].time = 1.0f;

			colourGradient.SetKeys(colourKeys, alphaKeys);
		}
	}

	void ResetSimulation()
	{
		// Calculate actual simulation size based on resolution multiplier
		currentSize = Mathf.RoundToInt(size * resolutionMultiplier);

		// Calculate cell size for physical accuracy
		cellSize = physicalSize / currentSize;

		// Scale time step inversely with resolution to maintain stability
		dtScale = autoAdjustParameters ? 128f / currentSize : 1f;

		// Initialize arrays
		int totalSize = currentSize * currentSize;
		density = new float[totalSize];
		velocityX = new float[totalSize];
		velocityY = new float[totalSize];
		velocityX0 = new float[totalSize];
		velocityY0 = new float[totalSize];

		// Create or recreate visualization texture
		if (fluidTexture != null)
		{
			Destroy(fluidTexture);
		}
		fluidTexture = new Texture2D(currentSize, currentSize);
		fluidTexture.filterMode = FilterMode.Bilinear;

		// Create or update material for rendering
		if (fluidMaterial == null)
		{
			fluidMaterial = new Material(Shader.Find("Unlit/Texture"));
		}
		fluidMaterial.mainTexture = fluidTexture;

		// Create quad for visualization if it doesn't exist
		GameObject quad = GameObject.Find("FluidSimulationQuad");
		if (quad == null)
		{
			quad = GameObject.CreatePrimitive(PrimitiveType.Quad);
			quad.name = "FluidSimulationQuad";
			quad.GetComponent<Renderer>().material = fluidMaterial;
		}
		else
		{
			quad.GetComponent<Renderer>().material = fluidMaterial;
		}

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

		Debug.Log($"Fluid simulation reset with resolution: {currentSize}x{currentSize}, cell size: {cellSize}");
	}

	void Update()
	{
		elapsedTime += Time.deltaTime;

		// Handle source positioning with mouse if enabled
		if (moveSourceWithMouse && Input.GetKey(sourcePositionKey))
		{
			Vector2 mousePos = GetMousePositionInGrid();
			sourcePositionX = Mathf.Clamp01(mousePos.x / currentSize);
			sourcePositionY = Mathf.Clamp01(mousePos.y / currentSize);
		}

		// Process custom source if enabled
		if (enableCustomSource)
		{
			UpdateCustomSource();
		}

		// Add forces based on mouse input (only when not positioning the source)
		if (Input.GetMouseButton(0) && !(moveSourceWithMouse && Input.GetKey(sourcePositionKey)))
		{
			Vector2 mousePos = GetMousePositionInGrid();
			if (mousePos.x >= 0 && mousePos.x < currentSize && mousePos.y >= 0 && mousePos.y < currentSize)
			{
				// Scale density and velocity with resolution
				float densityAmount = 100f * resolutionMultiplier;
				float velocityScale = 10f * resolutionMultiplier;

				AddDensity(mousePos.x, mousePos.y, densityAmount);
				AddVelocity(mousePos.x, mousePos.y,
					Input.GetAxis("Mouse X") * velocityScale,
					Input.GetAxis("Mouse Y") * velocityScale);
			}
		}

		Simulate();
		UpdateVisualization();
	}

	void UpdateCustomSource()
	{
		// Convert normalized position (0-1) to grid coordinates
		float sourceX = sourcePositionX * currentSize;
		float sourceY = sourcePositionY * currentSize;

		// Calculate effective strength for pulsing if enabled
		float effectiveStrength = sourceStrength;
		if (sourcePulsing)
		{
			// Create pulsing effect using sine wave
			float pulseScale = Mathf.Abs(Mathf.Sin(elapsedTime * sourcePulseRate * Mathf.PI));
			effectiveStrength *= pulseScale;
		}

		// Adjust strength based on resolution
		effectiveStrength *= resolutionMultiplier;

		// Apply density to source with radius
		float radiusInCells = sourceRadius * resolutionMultiplier;

		for (int i = Mathf.Max(0, Mathf.FloorToInt(sourceX - radiusInCells));
			 i <= Mathf.Min(currentSize - 1, Mathf.CeilToInt(sourceX + radiusInCells)); i++)
		{
			for (int j = Mathf.Max(0, Mathf.FloorToInt(sourceY - radiusInCells));
				 j <= Mathf.Min(currentSize - 1, Mathf.CeilToInt(sourceY + radiusInCells)); j++)
			{
				// Calculate distance from source point
				float distSq = (i - sourceX) * (i - sourceX) + (j - sourceY) * (j - sourceY);
				float dist = Mathf.Sqrt(distSq);

				// Only add density within radius
				if (dist <= radiusInCells)
				{
					// Reduce strength as we get farther from center (linear falloff)
					float falloff = 1.0f - (dist / radiusInCells);
					AddDensity(i, j, effectiveStrength * falloff);

					// Add velocity if enabled
					if (sourceEmitsVelocity)
					{
						// Convert direction angle to velocity components
						float angleRad = sourceDirection * Mathf.Deg2Rad;
						float vx = Mathf.Cos(angleRad) * sourceVelocity * resolutionMultiplier * falloff;
						float vy = Mathf.Sin(angleRad) * sourceVelocity * resolutionMultiplier * falloff;

						AddVelocity(i, j, vx, vy);
					}
				}
			}
		}
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

		return new Vector2(normalizedX * currentSize, normalizedY * currentSize);
	}

	void Simulate()
	{
		// Scale diffusion and viscosity with resolution
		float effectiveTimeStep = autoAdjustParameters ? timeStep * dtScale : timeStep;
		float effectiveDiffusion = autoAdjustParameters ? diffusion / resolutionMultiplier : diffusion;
		float effectiveViscosity = autoAdjustParameters ? viscosity / resolutionMultiplier : viscosity;

		VelocityStep(effectiveTimeStep, effectiveViscosity);
		DensityStep(effectiveTimeStep, effectiveDiffusion);
	}

	void VelocityStep(float dt, float visc)
	{
		Diffuse(1, velocityX0, velocityX, visc, dt);
		Diffuse(2, velocityY0, velocityY, visc, dt);

		Project(velocityX0, velocityY0, velocityX, velocityY);

		Advect(1, velocityX, velocityX0, velocityX0, velocityY0, dt);
		Advect(2, velocityY, velocityY0, velocityX0, velocityY0, dt);

		Project(velocityX, velocityY, velocityX0, velocityY0);
	}

	void DensityStep(float dt, float diff)
	{
		float[] densityTemp = new float[currentSize * currentSize];
		Diffuse(0, densityTemp, density, diff, dt);
		Advect(0, density, densityTemp, velocityX, velocityY, dt);
	}

	void AddDensity(float x, float y, float amount)
	{
		int i = Mathf.Clamp((int)x, 0, currentSize - 1);
		int j = Mathf.Clamp((int)y, 0, currentSize - 1);

		density[IX(i, j)] += amount;
	}

	void AddVelocity(float x, float y, float amountX, float amountY)
	{
		int i = Mathf.Clamp((int)x, 0, currentSize - 1);
		int j = Mathf.Clamp((int)y, 0, currentSize - 1);

		velocityX[IX(i, j)] += amountX;
		velocityY[IX(i, j)] += amountY;
	}

	void Diffuse(int b, float[] x, float[] x0, float diff, float dt)
	{
		float a = dt * diff * (currentSize - 2) * (currentSize - 2);
		LinearSolve(b, x, x0, a, 1 + 6 * a);
	}

	void LinearSolve(int b, float[] x, float[] x0, float a, float c)
	{
		for (int k = 0; k < 20; k++)
		{
			for (int i = 1; i < currentSize - 1; i++)
			{
				for (int j = 1; j < currentSize - 1; j++)
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
		for (int i = 1; i < currentSize - 1; i++)
		{
			for (int j = 1; j < currentSize - 1; j++)
			{
				div[IX(i, j)] = -0.5f * (
					velocX[IX(i + 1, j)] - velocX[IX(i - 1, j)] +
					velocY[IX(i, j + 1)] - velocY[IX(i, j - 1)]
				) / currentSize;
				p[IX(i, j)] = 0;
			}
		}

		SetBoundary(0, div);
		SetBoundary(0, p);
		LinearSolve(0, p, div, 1, 6);

		for (int i = 1; i < currentSize - 1; i++)
		{
			for (int j = 1; j < currentSize - 1; j++)
			{
				velocX[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * currentSize;
				velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * currentSize;
			}
		}

		SetBoundary(1, velocX);
		SetBoundary(2, velocY);
	}

	void Advect(int b, float[] d, float[] d0, float[] velocX, float[] velocY, float dt)
	{
		float dt0 = dt * (currentSize - 2);

		for (int i = 1; i < currentSize - 1; i++)
		{
			for (int j = 1; j < currentSize - 1; j++)
			{
				float x = i - dt0 * velocX[IX(i, j)];
				float y = j - dt0 * velocY[IX(i, j)];

				if (x < 0.5f) x = 0.5f;
				if (x > currentSize - 1.5f) x = currentSize - 1.5f;
				int i0 = (int)x;
				int i1 = i0 + 1;

				if (y < 0.5f) y = 0.5f;
				if (y > currentSize - 1.5f) y = currentSize - 1.5f;
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
		for (int i = 1; i < currentSize - 1; i++)
		{
			x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
			x[IX(currentSize - 1, i)] = b == 1 ? -x[IX(currentSize - 2, i)] : x[IX(currentSize - 2, i)];
			x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
			x[IX(i, currentSize - 1)] = b == 2 ? -x[IX(i, currentSize - 2)] : x[IX(i, currentSize - 2)];
		}

		x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
		x[IX(0, currentSize - 1)] = 0.5f * (x[IX(1, currentSize - 1)] + x[IX(0, currentSize - 2)]);
		x[IX(currentSize - 1, 0)] = 0.5f * (x[IX(currentSize - 2, 0)] + x[IX(currentSize - 1, 1)]);
		x[IX(currentSize - 1, currentSize - 1)] = 0.5f * (x[IX(currentSize - 2, currentSize - 1)] + x[IX(currentSize - 1, currentSize - 2)]);
	}

	int IX(int x, int y)
	{
		return x + y * currentSize;
	}

	void UpdateVisualization()
	{
		Color[] colours = new Color[currentSize * currentSize];

		for (int i = 0; i < currentSize; i++)
		{
			for (int j = 0; j < currentSize; j++)
			{
				float d = density[IX(i, j)] * colourIntensity;

				if (useGradient)
				{
					// Use gradient evaluation (d is clamped to 0-1 range)
					colours[IX(i, j)] = colourGradient.Evaluate(Mathf.Clamp01(d));
				}
				else
				{
					// Use single colour with varying intensity
					Color baseColour = fluidcolour;
					colours[IX(i, j)] = new Color(
						baseColour.r * d,
						baseColour.g * d,
						baseColour.b * d,
						baseColour.a
					);
				}

				// Visualize source position if enabled
				if (visualizeSourcePosition && enableCustomSource)
				{
					float sourceX = sourcePositionX * currentSize;
					float sourceY = sourcePositionY * currentSize;
					float visualMarkerRadius = 3; // Size of the position marker in pixels

					float distSq = (i - sourceX) * (i - sourceX) + (j - sourceY) * (j - sourceY);
					if (distSq < visualMarkerRadius * visualMarkerRadius)
					{
						// Mark the source position with a distinct color
						colours[IX(i, j)] = sourcePositionColor;
					}
				}
			}
		}

		fluidTexture.SetPixels(colours);
		fluidTexture.Apply();
	}

	// Helper method for debugging resolution
	public string GetCurrentResolution()
	{
		return $"{currentSize}x{currentSize} (Base: {size}, Multiplier: {resolutionMultiplier:F2})";
	}

	// Helper method to get the current source position	
	public Vector2 GetSourcePosition()
	{
		return new Vector2(sourcePositionX * currentSize, sourcePositionY * currentSize);
	}

	// Helper method to manually set the source position using grid coordinates
	public void SetSourcePosition(float x, float y)
	{
		sourcePositionX = Mathf.Clamp01(x / currentSize);
		sourcePositionY = Mathf.Clamp01(y / currentSize);
	}
}