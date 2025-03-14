using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Collections;
using Unity.Jobs;
using Unity.Burst;

namespace Assets.Scripts
{
    class Burst
    {
    }

	//[BurstCompile]
	public struct DiffuseJob : IJobParallelFor
	{
		[ReadOnly] public NativeArray<float> x0;
		public NativeArray<float> x;
		[ReadOnly] public NativeArray<bool> obstacles;

		public int size;
		public float a;
		public float c;
		public int b; // Boundary condition flag
		public int iteration; // Current iteration

		public void Execute(int index)
		{
			// Skip boundary cells (we'll handle those separately)
			int i = index % size;
			int j = index / size;

			if (i <= 0 || i >= size - 1 || j <= 0 || j >= size - 1)
				return;

			// Skip obstacles
			if (obstacles[index])
				return;

			// Same calculation as in the original LinearSolve
			x[index] = (x0[index] + a * (
				x[i + 1 + j * size] + x[i - 1 + j * size] +
				x[i + (j + 1) * size] + x[i + (j - 1) * size]
			)) / c;
		}
	}
}
