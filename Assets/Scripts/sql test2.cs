using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Data;
using UnityEngine;
using Mono.Data.Sqlite;

namespace Assets.Scripts
{
    class sql_test2 : MonoBehaviour
	{

		// Keep count of the clicks.
		[SerializeField] private int hitCount; // 1

		private void Start() // 2
		{
			// Read the persisted data and set the initial hit count.
			hitCount = 0; // 3
		}

		private void OnMouseDown() // 4
		{
			// Increment the hit count on each click and save the data.
			hitCount++; // 5
		}
	}
}
