using System.Collections;
using System.Collections.Generic;
using System.Data;
using UnityEngine;
using Mono.Data.Sqlite;
using static FluidSimulation;

public class sql_test : MonoBehaviour  
{
	// Start is called before the first frame update
	[SerializeField] private int hitCount = 0;

	void Start() // 13
	{
		// Read all values from the table.
		IDbConnection dbConnection = CreateAndOpenDatabase(); // 14
		IDbCommand dbCommandReadValues = dbConnection.CreateCommand(); // 15
		dbCommandReadValues.CommandText = "SELECT * FROM Chris"; // 16
		IDataReader dataReader = dbCommandReadValues.ExecuteReader(); // 17

		while (dataReader.Read())
		{
			// The `id` has index 0, our `hits` have the index 1.
			hitCount = dataReader.GetInt32(1); // 19
		}

		// Remember to always close the connection at the end.
		dbConnection.Close(); // 20
	}

	private void OnMouseDown()
	{
		hitCount++;

		// Insert hits into the table.
		IDbConnection dbConnection = CreateAndOpenDatabase(); // 2
		IDbCommand dbCommandInsertValue = dbConnection.CreateCommand(); // 9
		dbCommandInsertValue.CommandText = "INSERT OR REPLACE INTO PleaseWork (id, hits) VALUES (0, " + hitCount + ")"; // 10
		dbCommandInsertValue.ExecuteNonQuery(); // 11
		

		// Remember to always close the connection at the end.
		dbConnection.Close(); // 12
	}

	private IDbConnection CreateAndOpenDatabase() // 3
	{
		// Open a connection to the database.
		string dbUri = "URI=file:C:\\Users\\chris\\My project (2)\\test.db"; // 4
		IDbConnection dbConnection = new SqliteConnection(dbUri); // 5
		dbConnection.Open(); // 6

		// Create a table for the hit count in the database if it does not exist yet.
		IDbCommand dbCommandCreateTable = dbConnection.CreateCommand(); // 6
		dbCommandCreateTable.CommandText = "CREATE TABLE IF NOT EXISTS PleaseWork (id INTEGER PRIMARY KEY, hits INTEGER )"; // 7
		dbCommandCreateTable.ExecuteReader(); // 8

		return dbConnection;
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
