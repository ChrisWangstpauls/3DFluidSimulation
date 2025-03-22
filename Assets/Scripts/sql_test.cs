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

public class sql_test : MonoBehaviour
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



















	//public static void SaveSimulationData(int timeStep, FluidSimulation fluidSimulation)
	//{
	//	using var connection = new SqliteConnection("URI=file:C:\\Users\\Chris\\my project (2)\\test.db");
	//	connection.Open();

	//	for (int i = 0; i < fluidSimulation.currentSize; i++)
	//	{
	//		for (int j = 0; j < fluidSimulation.currentSize; j++)
	//		{
	//			using var command = connection.CreateCommand();
	//			int cellIndex = GridUtils.IX(i, j, fluidSimulation.currentSize);

	//			command.CommandText = "INSERT INTO SaveSimulationData(TimeStep, X, Y, Density, VelocityX, VelocityY, Pressure) " +
	//								  "VALUES (@timeStep, @x, @y, @density, @velocityX, @velocityY, @pressure)";

	//			command.Parameters.Clear();
	//			command.Parameters.AddWithValue("@timeStep", timeStep);
	//			command.Parameters.AddWithValue("@x", i);
	//			command.Parameters.AddWithValue("@y", j);
	//			command.Parameters.AddWithValue("@density", fluidSimulation.density[cellIndex]);
	//			command.Parameters.AddWithValue("@velocityX", fluidSimulation.velocityX[cellIndex]);
	//			command.Parameters.AddWithValue("@velocityY", fluidSimulation.velocityY[cellIndex]);
	//			command.Parameters.AddWithValue("@pressure", fluidSimulation.pressure[cellIndex]);

	//			command.ExecuteNonQuery();
	//		}
	//	}
	//	connection.Close();
	//}
}
