using System.Collections;
using System.Collections.Generic;
using System.Data;
using UnityEngine;
using Mono.Data.Sqlite;

using static FluidSimulation;
using System.Data.Common;
using System.Threading.Tasks;
using System;
using System.IO;

public class SQL : MonoBehaviour
{

	public static void SaveSimRunParams(int size, float diffusion, float viscosity, float timeStep,
	   bool sourceEnabled, float sourceStrength, float sourceX, float sourceY,
	   bool obstacleEnabled, string obstacleType, float obstacleX, float obstacleY,
	   float obstacleRadius, float obstacleWidth, float obstacleHeight)
	{
		using (var conn = new SqliteConnection("URI=file:C:\\Users\\chris\\My project (2)\\test.db"))
		{
			conn.Open();
			using (var cmd = conn.CreateCommand())
			{
				cmd.CommandText = @"INSERT INTO SimulationRuns 
                    (Size, Diffusion, Viscosity, TimeStep, SourceEnabled, SourceStrength, SourcePositionX, SourcePositionY, 
                     ObstacleEnabled, ObstacleType, ObstaclePositionX, ObstaclePositionY, ObstacleRadius, ObstacleWidth, ObstacleHeight)
                    VALUES 
                    (@Size, @Diffusion, @Viscosity, @TimeStep, @SourceEnabled, @SourceStrength, @SourcePositionX, @SourcePositionY, 
                     @ObstacleEnabled, @ObstacleType, @ObstaclePositionX, @ObstaclePositionY, @ObstacleRadius, @ObstacleWidth, @ObstacleHeight)";

				cmd.Parameters.Add(new SqliteParameter("@Size", size));
				cmd.Parameters.Add(new SqliteParameter("@Diffusion", diffusion));
				cmd.Parameters.Add(new SqliteParameter("@Viscosity", viscosity));
				cmd.Parameters.Add(new SqliteParameter("@TimeStep", timeStep));
				cmd.Parameters.Add(new SqliteParameter("@SourceEnabled", sourceEnabled ? 1 : 0));
				cmd.Parameters.Add(new SqliteParameter("@SourceStrength", sourceStrength));
				cmd.Parameters.Add(new SqliteParameter("@SourcePositionX", sourceX));
				cmd.Parameters.Add(new SqliteParameter("@SourcePositionY", sourceY));
				cmd.Parameters.Add(new SqliteParameter("@ObstacleEnabled", obstacleEnabled ? 1 : 0));
				cmd.Parameters.Add(new SqliteParameter("@ObstacleType", obstacleType));
				cmd.Parameters.Add(new SqliteParameter("@ObstaclePositionX", obstacleX));
				cmd.Parameters.Add(new SqliteParameter("@ObstaclePositionY", obstacleY));
				cmd.Parameters.Add(new SqliteParameter("@ObstacleRadius", obstacleRadius));
				cmd.Parameters.Add(new SqliteParameter("@ObstacleWidth", obstacleWidth));
				cmd.Parameters.Add(new SqliteParameter("@ObstacleHeight", obstacleHeight));

				cmd.ExecuteNonQuery();

			}
			conn.Close();
		}
	}

	public static void LogRuntimeMetrics(
		int step,
		float avgDensity,
		float maxVelocity,
		float sourceX,
		float sourceY,
		float sourceStrength,
		bool obstacleEnabled)
	{
		using (var conn = new SqliteConnection("URI=file:C:\\Users\\chris\\My project (2)\\test.db"))
		{
			conn.Open();
			using (var cmd = conn.CreateCommand())
			{
				cmd.CommandText = @"INSERT INTO RuntimeMetrics 
                (Step, AverageDensity, MaxVelocityMagnitude, CurrentSourceX, 
                 CurrentSourceY, CurrentSourceStrength, ObstaclePresent)
                VALUES 
                (@step, @avgDensity, @maxVelocity, @sourceX, 
                 @sourceY, @sourceStrength, @obstacle)";

				cmd.Parameters.AddRange(new[] {
				new SqliteParameter("@step", step),
				new SqliteParameter("@avgDensity", avgDensity),
				new SqliteParameter("@maxVelocity", maxVelocity),
				new SqliteParameter("@sourceX", sourceX),
				new SqliteParameter("@sourceY", sourceY),
				new SqliteParameter("@sourceStrength", sourceStrength),
				new SqliteParameter("@obstacle", obstacleEnabled ? 1 : 0)
			});

				cmd.ExecuteNonQuery();
			}
			conn.Close();
		}
	}
}


