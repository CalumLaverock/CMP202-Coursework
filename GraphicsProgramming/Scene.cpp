// Data Structures and Algorithms II : Intro to AMP and benchmarking exercise
// Ruth Falconer  <r.falconer@abertay.ac.uk>
// Adapted from C++ AMP book http://ampbook.codeplex.com/license.

#include "Scene.h"

#define WIDTH 1280
#define HEIGHT 960
#define MAX_ITERATIONS 1080

// Need to access the concurrency libraries 
using namespace concurrency;

// Import things we need from the standard library
using std::chrono::duration_cast;
using std::chrono::microseconds;
using std::cout;
using std::cout;
using std::endl;

// Define the alias "the_clock" for the clock type we're going to use.
typedef std::chrono::steady_clock the_serial_clock;
typedef std::chrono::steady_clock the_amp_clock;

uint32_t image[HEIGHT][WIDTH];

Scene::Scene(Input* in)
{
	input = in;

	std::vector<concurrency::accelerator> accls = accelerator::get_all();

	for (concurrency::accelerator a : accls)
	{
		if (!a.is_emulated)
		{
			usable_accls.push_back(a);
		}
	}

	accelerator::set_default(usable_accls[0].device_path);

	// Check AMP support
	//query_AMP_support();

	test_left = -2; 
	test_right = 1;
	test_top = -1.25;
	test_bottom = 1.25;

	RunMandelbrot((-2.5 * zoom) + pan_x, (1.5 * zoom) + pan_x, (-2.5 * zoom) + pan_y, (1.5 * zoom) + pan_y);
}

void Scene::report_accelerator(const accelerator a)
{
	const std::wstring bs[2] = { L"false", L"true" };
	std::wcout << ": " << a.description << " "
		<< endl << "       device_path                       = " << a.device_path
		<< endl << "       dedicated_memory                  = " << std::setprecision(4) << float(a.dedicated_memory) / (1024.0f * 1024.0f) << " Mb"
		<< endl << "       has_display                       = " << bs[a.has_display]
		<< endl << "       is_debug                          = " << bs[a.is_debug]
		<< endl << "       is_emulated                       = " << bs[a.is_emulated]
		<< endl << "       supports_double_precision         = " << bs[a.supports_double_precision]
		<< endl << "       supports_limited_double_precision = " << bs[a.supports_limited_double_precision]
		<< endl;
}

// List and select the accelerator to use
void Scene::list_accelerators()
{
	//get all accelerators available to us and store in a vector so we can extract details
	std::vector<accelerator> accls = accelerator::get_all();

	// iterates over all accelerators and print characteristics
	for (int i = 0; i < accls.size(); i++)
	{
		accelerator a = accls[i];
		report_accelerator(a);

	}

	accelerator acc = accelerator(accelerator::default_accelerator);
	std::wcout << " default acc = " << acc.description << endl;

} // list_accelerators

  // query if AMP accelerator exists on hardware
void Scene::query_AMP_support()
{
	std::vector<accelerator> accls = accelerator::get_all();
	if (accls.empty())
	{
		cout << "No accelerators found that are compatible with C++ AMP" << std::endl;
	}
	else
	{
		cout << "Accelerators found that are compatible with C++ AMP" << std::endl;
		list_accelerators();
	}
} // query_AMP_support

void Scene::TaskDecomposition(int num_of_tasks, double left_, double right_, double top_, double bottom_)
{
	int index = 0;
	rows_per_task_ = HEIGHT / num_of_tasks;
	int start = 0;

	double top_inc = abs((top_ - bottom_) / num_of_tasks);
	double top_temp = top_;
	double bottom_temp = top_ + top_inc;

	for (int i = 0; i < num_of_tasks; i++)
	{
		TaskData* t = new TaskData(usable_accls[0]);
		t->id = index;
		index++;
		t->start_row = start;
		t->write_ext = { rows_per_task_, WIDTH };
		t->left = left_;
		t->right = right_;
		t->top = top_temp;
		t->bottom = bottom_temp;

		start += rows_per_task_;
		top_temp += top_inc;
		bottom_temp = top_temp + top_inc;

		tasks.push(t);
	}
}

//HSV implementation used from Peter Hall git repository given to use for
//CMP208: https://gitlab.com/cmp208-2018/polishing_a_game/commit/6d3a58cb32497df257fc6b3bb6381b8c9c65d5d4
//both functions, "makeRGB" and "colourHSVtoRGB" from that repository
uint32_t makeRGB(double r, double g, double b) restrict(amp)
{
	return (int)(r * 255) | (int)(g * 255) << 8 | (int)(b * 255) << 16;
}

uint32_t colourHSVtoRGB(double h, double s, double v) restrict(amp)
{
	double i;
	while (h >= 360.0f)
		h -= 360.0f;
	h /= 60.f;

	//need to use (i = (double)floor(h)) but can't because of restrict(amp)
	//so using this to round h down by casting to int to cut off decimal point
	//then casting back to double as rounded down value.
	int h_int = (int)h;
	i = (double)h_int;

	//need to subtract 1 after cutting off decimal value to round negative numbers down.
	if (i < 0)
	{
		i -= 1.f;
	}

	double f = h - i;
	double p = v * (1 - s);
	double q = v * (1 - (s * f));
	double t = v * (1 - (s * (1 - f)));

	switch ((int)i)
	{
	case 0: return makeRGB(v, t, p);
	case 1: return makeRGB(q, v, p);
	case 2: return makeRGB(p, v, t);
	case 3: return makeRGB(p, q, v);
	case 4: return makeRGB(t, p, v);
	case 5: return makeRGB(v, p, q);
	}
	return 0;
}

struct Complex
{
	double r, i;
};

Complex c_add(Complex c1, Complex c2) restrict(amp)
{
	return
	{
		//real part
		c1.r + c2.r,
		//imaginary part
		c1.i + c2.i
	};
}

Complex c_mult(Complex c1, Complex c2) restrict(amp)
{
	return
	{
		//real part
		((c1.r * c2.r) - (c1.i * c2.i)),
		//imaginary part
		((c1.r * c2.i) + (c1.i * c2.r))
	};
}

double c_abs(Complex c) restrict(amp)
{
	return (c.r * c.r) + (c.i * c.i);
}

void Scene::WriteTGA(const char *filename)
{
	std::ofstream outfile(filename, std::ofstream::binary);

	uint8_t header[18] = {
		0, // no image ID
		0, // no colour map
		2, // uncompressed 24-bit image
		0, 0, 0, 0, 0, // empty colour map specification
		0, 0, // X origin
		0, 0, // Y origin
		WIDTH & 0xFF, (WIDTH >> 8) & 0xFF, // width
		HEIGHT & 0xFF, (HEIGHT >> 8) & 0xFF, // height
		24, // bits per pixel
		0, // image descriptor
	};
	outfile.write((const char *)header, 18);

	for (int y = 0; y < HEIGHT; ++y)
	{
		for (int x = 0; x < WIDTH; ++x)
		{
			uint8_t pixel[3] = {
				(image[y][x] >> 16) & 0xFF, // blue channel
				(image[y][x] >> 8) & 0xFF, // green channel
				(image[y][x]) & 0xFF, // red channel
			};
			outfile.write((const char *)pixel, 3);
		}
	}

	outfile.close();
	if (!outfile)
	{
		// An error has occurred at some point since we opened the file.
		std::cout << "Error writing to " << filename << std::endl;
		exit(1);
	}
}

std::chrono::microseconds Scene::RunMandelbrot(double left_, double right_, double top_, double bottom_, int num_tasks)
{
	TaskDecomposition(num_tasks, left_, right_, top_, bottom_);

	accl_queue.clear();

	for (concurrency::accelerator a : usable_accls)
	{
		accl_queue.push(a);
	}

	TaskData* task;
	concurrency::accelerator accl;

	start = std::chrono::high_resolution_clock::now();

	//compute mandelbrot using multiple accelerators if testing with multi GPU
	//otherwise only run on single GPU
	if (multi_testing)
	{
		//if testing on multi GPU, pop a task and accelerator off of their respective
		//concurrent queues and pass them into ComputeMandelbrot
		while (tasks.try_pop(task) && accl_queue.try_pop(accl))
		{
			ComputeMandelbrot(task, accl);
			accl_queue.push(accl);
		}
	}
	else
	{
		//if running on single GPU, pop a task off of the accelerator
		while (tasks.try_pop(task))
		{
			ComputeMandelbrot(task);
		}
	}

	for (int i = 0; i < threads.size(); i++)
	{
		threads[i]->join();
		delete threads[i];
	}

	threads.clear();

	end = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<microseconds>(end - start);

	return duration;
}

void Scene::ComputeMandelbrot(TaskData* t)
{
	threads.push_back(new std::thread([&]() {
		double left = t->left;
		double right = t->right;
		double top = t->top;
		double bottom = t->bottom;

		int rows_per_task = rows_per_task_;

		concurrency::array_view<uint32_t, 2> t_image(t->write_ext, (uint32_t*)image[t->start_row]);
		t_image.discard_data();

		// It is wise to use exception handling here - AMP can fail for many reasons
		// and it useful to know why (e.g. using double precision when there is limited or no support)
		try
		{
			concurrency::parallel_for_each(t->view, t->write_ext, [=](concurrency::index<2> idx)  restrict(amp)
			{
				int x = idx[1];
				int y = idx[0];

				// Work out the point in the complex plane that
				// corresponds to this pixel in the output image.
				Complex c =
				{
					left + (x * (right - left) / WIDTH),
					top + (y * (bottom - top) / rows_per_task)
				};

				// Start off z at (0, 0).
				Complex z =
				{
					0.0,
					0.0
				};

				// Iterate z = z^2 + c until z moves more than 2 units
				// away from (0, 0), or we've iterated too many times.
				int iterations = 0;
				//checking if within bounds of 2^2 as c_abs doesn't square root for perfomance reasons
				//so the complex number and the bounds are squared meaning outcome is the same
				while (c_abs(z) < 4.0 && iterations < MAX_ITERATIONS)
				{
					z = c_add(c_mult(z, z), c);

					++iterations;
				}

				if (iterations == MAX_ITERATIONS)
				{
					// z didn't escape from the circle.
					// This point is in the Mandelbrot set.
					t_image[idx] = 0x000000; // black
				}
				else
				{
					// z escaped within less than MAX_ITERATIONS
					// iterations. This point isn't in the set.
					t_image[idx] = colourHSVtoRGB(iterations, 0.8f, 0.8f); // SHIFT
				}
			});
			t_image.synchronize();
		}
		catch (const Concurrency::runtime_exception& ex)
		{
			MessageBoxA(NULL, ex.what(), "Error", MB_ICONERROR);
		}
	}));

	delete t;
}

void Scene::ComputeMandelbrot(TaskData* t, concurrency::accelerator a)
{
	threads.push_back(new std::thread([&]() {
		double left = t->left;
		double right = t->right;
		double top = t->top;
		double bottom = t->bottom;

		t->view = a.default_view;

		//std::wcout << a.description << std::endl;

		int rows_per_task = rows_per_task_;

		concurrency::array_view<uint32_t, 2> t_image(t->write_ext, (uint32_t*)image[t->start_row]);
		t_image.discard_data();

		// It is wise to use exception handling here - AMP can fail for many reasons
		// and it useful to know why (e.g. using double precision when there is limited or no support)
		try
		{
  			concurrency::parallel_for_each(t->view, t->write_ext, [=](concurrency::index<2> idx)  restrict(amp)
			{
				int x = idx[1];
				int y = idx[0];

				// Work out the point in the complex plane that
				// corresponds to this pixel in the output image.
				Complex c =
				{
					left + (x * (right - left) / WIDTH),
					top + (y * (bottom - top) / rows_per_task)
				};

				// Start off z at (0, 0).
				Complex z =
				{
					0.0,
					0.0
				};

				// Iterate z = z^2 + c until z moves more than 2 units
				// away from (0, 0), or we've iterated too many times.
				int iterations = 0;
				//checking if within bounds of 2^2 as c_abs doesn't square root for perfomance reasons
				//so the complex number and the bounds are squared meaning outcome is the same
				while (c_abs(z) < 4.0 && iterations < MAX_ITERATIONS)
				{
					z = c_add(c_mult(z, z), c);

					++iterations;
				}

				if (iterations == MAX_ITERATIONS)
				{
					// z didn't escape from the circle.
					// This point is in the Mandelbrot set.
					t_image[idx] = 0x000000; // black
				}
				else
				{
					// z escaped within less than MAX_ITERATIONS
					// iterations. This point isn't in the set.
					t_image[idx] = colourHSVtoRGB(iterations, 0.8f, 0.8f); // SHIFT
				}
			});
			t_image.synchronize();
		}
		catch (const Concurrency::runtime_exception& ex)
		{
			MessageBoxA(NULL, ex.what(), "Error", MB_ICONERROR);
		}
	}));

	delete t;
}

void Scene::update()
{
	if (!testing)
	{
		//toggle between multi GPU and single GPU
		if (input->isKeyDown('m'))
		{
			if (multi == false)
			{
				multi = true;
				std::cout << "Single GPU Testing Disabled\nMulti GPU Testing Enabled" << std::endl;
			}
			else
			{
				multi = false;
				std::cout << "Single GPU Testing Enabled\nMulti GPU Testing Disabled" << std::endl;
			}

			input->SetKeyUp('m');
		}

		//change which portion of mandelbrot to test to the current image
		if (input->isKeyDown('z'))
		{
			/*test_left = (-2.5 * zoom) + pan_x;
			test_right = (1.5 * zoom) + pan_x;
			test_top = (-2.5 * zoom) + pan_y;
			test_bottom = (1.5 * zoom) + pan_y;*/

			std::cout << "New Testing Area: " << std::endl;
			std::cout << test_left << ", " << test_right << ", " << test_top << ", " << test_bottom << ", " << std::endl;
			input->SetKeyUp('z');
		}

		//'w', 'a', 's' and 'd' to pan around Mandelbrot
		if (input->isKeyDown('a'))
		{
			input->SetKeyUp('a');

			//multiply by (zoom * zoom_increment / 2) in order to scale the panning amount as the
			//mandelbrot zooms in. Without this, panning while zoomed jumps by large amounts
			pan_x -= ((0.5 / 120.f) * (zoom * zoom_increment / 2));

			RunMandelbrot((-2.5 * zoom) + pan_x, (1.5 * zoom) + pan_x, (-2.5 * zoom) + pan_y, (1.5 * zoom) + pan_y);
		}
		else if (input->isKeyDown('d'))
		{
			input->SetKeyUp('d');

			pan_x += ((0.5 / 120.f) * (zoom * zoom_increment / 2));
			
			RunMandelbrot((-2.5 * zoom) + pan_x, (1.5 * zoom) + pan_x, (-2.5 * zoom) + pan_y, (1.5 * zoom) + pan_y);
		}
		else if (input->isKeyDown('w'))
		{
			input->SetKeyUp('w');

			pan_y += ((0.5 / 120.f) * (zoom * zoom_increment / 2));
			
			RunMandelbrot((-2.5 * zoom) + pan_x, (1.5 * zoom) + pan_x, (-2.5 * zoom) + pan_y, (1.5 * zoom) + pan_y);
		}
		else if (input->isKeyDown('s'))
		{
			input->SetKeyUp('s');

			pan_y -= ((0.5 / 120.f) * (zoom * zoom_increment / 2));
			
			RunMandelbrot((-2.5 * zoom) + pan_x, (1.5 * zoom) + pan_x, (-2.5 * zoom) + pan_y, (1.5 * zoom) + pan_y);
		}

		//spacebar to zoom in, 'b' to zoom out
		if (input->isKeyDown(' '))
		{
			input->SetKeyUp(' ');

			zoom /= 1.02;
			zoom_increment++;
			
			RunMandelbrot((-2.5 * zoom) + pan_x, (1.5 * zoom) + pan_x, (-2.5 * zoom) + pan_y, (1.5 * zoom) + pan_y);
		}

		if (input->isKeyDown('b'))
		{
			input->SetKeyUp('b');

			zoom *= 1.02;
			zoom_increment--;
			
			RunMandelbrot((-2.5 * zoom) + pan_x, (1.5 * zoom) + pan_x, (-2.5 * zoom) + pan_y, (1.5 * zoom) + pan_y);
		}

		//'r' will reset the Mandelbrot back to it's original position
		if (input->isKeyDown('r'))
		{
			input->SetKeyUp('r');

			pan_x = 0.5, pan_y = 0.5;
			zoom = 1;
			zoom_increment = 1;

			RunMandelbrot((-2.5 * zoom) + pan_x, (1.5 * zoom) + pan_x, (-2.5 * zoom) + pan_y, (1.5 * zoom) + pan_y);
		}

		//'c' will take a screenshot and write it to a tga file
		if (input->isKeyDown('c'))
		{
			input->SetKeyUp('c');

			WriteTGA("Mandelbrot.tga");
		}

		//'t' will begin testing for multi GPU or single GPU depending on which is selected
		//at the time of 't' being pressed
		if (input->isKeyDown('t'))
		{
			input->SetKeyUp('t');

			multi_testing = multi ? true : false;

			testing = true;
		}
	}
	else
	{
		std::cout << "Testing " << (multi_testing ? "Multi GPU" : "Single GPU") << std::endl;

		PerformanceTestFarm();

		std:cout << "Stopped testing" << std::endl;
	}
}

void Scene::render()
{
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &texture);

	glBindTexture(GL_TEXTURE_2D, texture);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, WIDTH, HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);

	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glBegin(GL_QUADS);
	glTexCoord2i(0, 0);
	glVertex2f(-1.f, -1.f);

	glTexCoord2i(0, 1);
	glVertex2f(-1.f, 1.f);

	glTexCoord2i(1, 1);
	glVertex2f(1.f, 1.f);

	glTexCoord2i(1, 0);
	glVertex2f(1.f, -1.f);
	glEnd();

	glutSwapBuffers();
}

void Scene::PerformanceTestFarm()
{
	std::ofstream file;
	file.open("MandelbrotPerformanceTest.csv");

	for (int num_tasks = 2; num_tasks < 32; num_tasks *= 2)
	{
		std::vector<std::chrono::microseconds> times;

		for (iterations = 0; iterations < 51; iterations++)
		{
			//times.push_back(RunMandelbrot(test_left , test_right, test_top, test_bottom, num_tasks));
			times.push_back(RunMandelbrot(0.306234641146628976, 0.35456521296873689, 0.50062198797240853, 0.54895255947485566, num_tasks));
		}

		times.erase(times.begin());

		all_times.push_back(times);
	}

	for (int i = 0; i < all_times.size(); i++)
	{
		for (int j = 0; j < all_times[i].size(); j++)
		{
			file << all_times[i][j].count() << std::endl;

			for (int k = 0; k < i; k++)
			{
				file << ",";
			}
		}
	}

	file.close();

	testing = false;
	multi_testing = false;

	RunMandelbrot((-2.5 * zoom) + pan_x, (1.5 * zoom) + pan_x, (-2.5 * zoom) + pan_y, (1.5 * zoom) + pan_y);
}