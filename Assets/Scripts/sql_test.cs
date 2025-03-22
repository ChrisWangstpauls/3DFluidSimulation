using System.Collections;
using System.Collections.Generic;
using System.Data;
using UnityEngine;
using Mono.Data.Sqlite;

using static FluidSimulation;
using System.Data.Common;
using System.Threading.Tasks;

public class sql_test : FluidSimulation
{
	public static void SaveSimulationData(int timeStep, FluidSimulation fluidSimulation)
	{
		using var connection = new SqliteConnection("URI=file:C:\\Users\\Chris\\my project (2)\\test.db");
		connection.Open();
		using var command = connection.CreateCommand();
		command.CommandText = "INSERT INTO SaveSimulationData(TimeStep) " +
					  "VALUES (@timeStep)";

		command.Parameters.Clear();
		command.Parameters.AddWithValue("@timeStep", timeStep);

		//for (int i = 0; i < fluidSimulation.currentSize; i++)
		//{
		//	for (int j = 0; j < fluidSimulation.currentSize; j++)
		//	{
		//		using var command = connection.CreateCommand();
		//		int cellIndex = GridUtils.IX(i, j, fluidSimulation.currentSize);

		//		command.CommandText = "INSERT INTO SaveSimulationData(TimeStep, X, Y, Density, VelocityX, VelocityY, Pressure) " +
		//							  "VALUES (@timeStep, @x, @y, @density, @velocityX, @velocityY, @pressure)";

		//		command.Parameters.Clear();
		//		command.Parameters.AddWithValue("@timeStep", timeStep);
		//		command.Parameters.AddWithValue("@x", i);
		//		command.Parameters.AddWithValue("@y", j);
		//		command.Parameters.AddWithValue("@density", fluidSimulation.density[cellIndex]);
		//		command.Parameters.AddWithValue("@velocityX", fluidSimulation.velocityX[cellIndex]);
		//		command.Parameters.AddWithValue("@velocityY", fluidSimulation.velocityY[cellIndex]);
		//		command.Parameters.AddWithValue("@pressure", fluidSimulation.pressure[cellIndex]);

		//		command.ExecuteNonQuery();
		//	}
		//}
		connection.Close();
	}
}
