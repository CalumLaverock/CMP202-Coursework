// Scene class. Configures a basic 3D scene.
// Interfaces with the Input class to handle user input
// Calculates and outputs Frames Per Second (FPS) rendered.
// Important functions are the constructor (initialising the scene), 
// update (for process user input and updating scene objects) and render (renders scene).
#ifndef _SCENE_H
#define _SCENE_H

#include <glut.h>
#include <gl/gl.h>
#include <gl/glu.h>

#include <chrono>
#include <iostream>
#include <iomanip>
#include <amp.h>
#include <time.h>
#include <string>
#include <array>
#include <amp_math.h>
#include <fstream>
#include <queue>
#include <concurrent_queue.h>
#include <thread>
#include <mutex>

#include "Input.h"

class Scene{

public:
	Scene(Input*);

	struct TaskData
	{
		int id;
		concurrency::accelerator_view view;
		int start_row;
		concurrency::extent<2> write_ext;
		double left, right, top, bottom;

		TaskData(concurrency::accelerator a) : view(a.default_view) {};
	};

	void render();
	void update();

private:
	void report_accelerator(const Concurrency::accelerator a);
	void list_accelerators();
	void query_AMP_support();

	void TaskDecomposition(int num_of_tasks, double left, double right, double top, double bottom);

	void WriteTGA(const char *filename);
	std::chrono::microseconds RunMandelbrot(double left, double right, double top, double bottom, int num_tasks = 2);

	void ComputeMandelbrot(TaskData*);
	void ComputeMandelbrot(TaskData*, concurrency::accelerator);
	void PerformanceTestFarm();

	Input* input;

	float pan_x = 0.5, pan_y = 0.5;
	float zoom = 1;
	float zoom_increment = 1;

	int rows_per_task_;

	std::vector<concurrency::accelerator> usable_accls;
	concurrency::concurrent_queue<concurrency::accelerator> accl_queue;
	concurrency::concurrent_queue<TaskData*> tasks;
	std::vector<std::thread*> threads;

	std::mutex task_mutex;

	//testing variables
	std::vector<std::vector<std::chrono::microseconds>> all_times;
	bool testing = false, multi = false, multi_testing = false;
	int iterations = 0;

	double test_left, test_right, test_top, test_bottom;
	bool test_zoomed = false;

	std::chrono::time_point<std::chrono::high_resolution_clock> start, end;

	//OpenGL variables
	GLuint texture;
};

#endif