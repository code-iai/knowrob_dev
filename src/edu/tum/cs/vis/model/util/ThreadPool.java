/*******************************************************************************
 * Copyright (c) 2012 Stefan Profanter. All rights reserved. This program and the accompanying
 * materials are made available under the terms of the GNU Public License v3.0 which accompanies
 * this distribution, and is available at http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors: Stefan Profanter - initial API and implementation, Year: 2012
 ******************************************************************************/
package edu.tum.cs.vis.model.util;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Simple Thread execution pool to execute multiple callable objects at same time.
 * 
 * @author Stefan Profanter
 * 
 */
public class ThreadPool {
	/**
	 * Main pool
	 */
	public static ExecutorService	pool;

	/**
	 * Executes the given callable objects in a thread pool and returns when all threads have
	 * finished and all callable objects have been executed.
	 * 
	 * @param threads
	 *            list of callable objects
	 */
	public static void executeInPool(List<Callable<Void>> threads) {

		if (pool == null)
			init();

		try {
			List<Future<Void>> futures = pool.invokeAll(threads);
			for (Future<Void> f : futures) {
				try {
					f.get(); // If called thread threw an exception, get will throw an
								// ExecutionException
				} catch (ExecutionException ex) {
					ex.getCause().printStackTrace(); // Print exceptions from called threads if
														// there were any
				}
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		threads.clear();
	}

	/**
	 * Initializes ExecutorService to use maximum ProcessorCount*25 threads at same time
	 */
	public static void init() {
		int threadNum = Runtime.getRuntime().availableProcessors() * 25;
		pool = Executors.newFixedThreadPool(threadNum);
	}
}
