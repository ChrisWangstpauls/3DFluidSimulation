using System.Collections;
using System.Collections.Generic;
using System.Data;
using UnityEngine;
using Mono.Data.Sqlite;
using static FluidSimulation;

public class sql_test : FluidSimulation  
{
    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        
    }

	/*void SaveSimulationData(int runID, int timeStep)
	{
		using (var connection = new SqliteConnection("URI=file:fluid_simulation.db"))
		{
			connection.Open();
			using (var command = connection.CreateCommand())
			{
				for (int i = 0; i < currentSize; i++)
				{
					for (int j = 0; j < currentSize; j++)
					{
						int cellIndex = GridUtils.IX(i, j, currentSize);

						command.CommandText = "INSERT INTO SimulationData (RunID, TimeStep, CellX, CellY, Density, VelocityX, VelocityY, Pressure) " +
											  "VALUES (@runID, @timeStep, @x, @y, @density, @velocityX, @velocityY, @pressure)";

						command.Parameters.Clear();
						command.Parameters.AddWithValue("@runID", runID);
						command.Parameters.AddWithValue("@timeStep", timeStep);
						command.Parameters.AddWithValue("@x", i);
						command.Parameters.AddWithValue("@y", j);
						command.Parameters.AddWithValue("@density", density[cellIndex]);
						command.Parameters.AddWithValue("@velocityX", velocityX[cellIndex]);
						command.Parameters.AddWithValue("@velocityY", velocityY[cellIndex]);
						command.Parameters.AddWithValue("@pressure", pressure[cellIndex]);

						command.ExecuteNonQuery();
					}
				}
			}
		}
	}*/

}
