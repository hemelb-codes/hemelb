// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef HEMELB_CFG_ON_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif// HEMELB_CFG_OSX
#include <rpc/types.h>
#include <rpc/rpc.h>
#include <cstring>
#include <string>

//  g++ -framework OpenGL -framework GLUT -framework Foundation visualize_images.c -o visualize_images

#define RAINBOW        0
#define HALF_RAINBOW   1
#define GREY           2

char *input_path;
char *output_path;

unsigned char *image_data, *image_buffer;

int pixels_x, pixels_y;
int images;
int image_count = 0;
int cycle_id = 1;

int min(int a, int b) {
	if (a < b) {
		return a;
	} else {
		return b;
	}
}

int max(int a, int b) {
	if (a > b) {
		return a;
	} else {
		return b;
	}
}

void visRainbowPalette(int id, float col[3]) {
	if (id == 0) {
		col[0] = 0.F;
		col[1] = 0.F;
		col[2] = 1.F;
	} else if (id == 1) {
		col[0] = 0.F;
		col[1] = 1.F;
		col[2] = 1.F;
	} else if (id == 2) {
		col[0] = 0.F;
		col[1] = 1.F;
		col[2] = 0.F;
	} else if (id == 3) {
		col[0] = 1.F;
		col[1] = 1.F;
		col[2] = 0.F;
	} else if (id == 4) {
		col[0] = 1.F;
		col[1] = 0.F;
		col[2] = 0.F;
	}
}

void OpenWindow(int pixels_x, int pixels_y) {
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(pixels_x - 1, pixels_y - 1);

	glutCreateWindow(" ");

	glDisable (GL_DEPTH_TEST);
	glDisable (GL_BLEND);
	glShadeModel (GL_SMOOTH);
	glDisable (GL_DITHER);

	glClear (GL_COLOR_BUFFER_BIT);
}

void Projection(int pixels_x, int pixels_y) {
	glLoadIdentity();
	gluOrtho2D(0.F, (float) (pixels_x - 1),
	0.F, (float)(pixels_y - 1));

glClearColor	(1.0F, 1.0F, 1.0F, 0.0F);
}

void DisplayFrame(int pixels_x, int pixels_y) {
	for (int j = 0; j < pixels_y; j++) {
		for (int i = 0; i < pixels_x; i++) {
			if (i > 3 && i <= pixels_x - 3 && j > 3 && j <= pixels_y - 3)
				continue;

			image_data[(j * pixels_x + i) * 3 + 0] = 0;
			image_data[(j * pixels_x + i) * 3 + 1] = 0;
			image_data[(j * pixels_x + i) * 3 + 2] = 0;
		}
	}
}

void visVisualiseString(float r, float g, float b, float x, float y,
		char *string, void *font) {
	glColor3f(r, g, b);
	glRasterPos2f(x, y);

	for (int i = 0; i < (int) strlen(string); i++) {
		glutBitmapCharacter(font, string[i]);
	}
	glEnd();
}

void visVisualiseBarPars(char *unit, float threshold_min, float threshold_max,
		float x_min, float y_min, float x_max) {
	float dx = 0.25F * (x_max - x_min);

	char string[256];

	for (int i = 0; i <= 4; i++) {
		sprintf(string, "%.2f",
				threshold_min + 0.25 * i * (threshold_max - threshold_min));

		visVisualiseString(0.F, 0.F, 0.F, x_min + i * dx - 20.0F, y_min - 16.0F,
				string, GLUT_BITMAP_HELVETICA_18);
	}
	visVisualiseString(0.F, 0.F, 0.F, x_max + 45, y_min - 16.0F, unit,
			GLUT_BITMAP_HELVETICA_18);
}

void visVisualiseColourBar(int type, float x_min, float y_min, float x_max,
		float y_max) {
	float dx = 0.25F * (x_max - x_min);
	float col[3], grey;

	int edge = 0;

	glBegin (GL_QUAD_STRIP);

	for (int i = 0; i <= 4; i++) {
		if (type == RAINBOW)
		{
			visRainbowPalette(edge++, col);
			glColor3fv(col);
		} else if (type == HALF_RAINBOW)
		{
			visRainbowPalette(edge++, col);

			for (int l = 0; l < 3; l++)
				col[l] = 0.5F * (1.0F + col[l]);
			glColor3fv(col);
		} else {
			grey = edge++ / (float) 8;
			glColor3f(grey, grey, grey);
		}
		glVertex2f(x_min + i * dx, y_min);
		glVertex2f(x_min + i * dx, y_max);
	}
	glEnd();

	glBegin (GL_LINES);

	for (int i = 0; i <= 4; i++) {
		if (type == GREY)
		{
			glColor3f(1.0F, 1.0F, 1.0F);
		} else {
			glColor3f(0.0F, 0.0F, 0.0F);
		}
		glVertex2f(x_min + i * dx, y_min);
		glVertex2f(x_min + i * dx, y_max);
	}
	glEnd();
}

void SaveWindowImage(int pixels_x, int pixels_y, char *file_name) {
	FILE *ppm_image_file_ptr = fopen(file_name, "wb");

	int i, j;

	unsigned char *data_p = NULL;

	glReadBuffer (GL_FRONT);

	data_p = image_data;

	for (j = 0; j < pixels_y; j++) {
		glReadPixels(0, j, pixels_x, 1, GL_RGB, GL_UNSIGNED_BYTE, image_buffer);

		for (i = 0; i < pixels_x; i++) {
			*data_p = image_buffer[i * 3];
			data_p++;
			*data_p = image_buffer[i * 3 + 1];
			data_p++;
			*data_p = image_buffer[i * 3 + 2];
			data_p++;
		}
	}
	fprintf(ppm_image_file_ptr, "P6\n%i %i\n255\n", pixels_x, pixels_y);

	for (j = pixels_y - 1; j >= 0; j--) {
		fwrite(image_data + j * pixels_x * 3, 1, pixels_x * 3,
				ppm_image_file_ptr);
	}
	fclose(ppm_image_file_ptr);
}

void Display(void) {
	FILE *input_image_file;
	FILE *temp_file_ptr;
	XDR xdr_image_file;

	float pressure_threshold_min, pressure_threshold_max;
	float velocity_threshold_max, stress_threshold_max;

	int mode;
	int bits_per_char = sizeof(char) * 8;
	int colour_mask = (1 << (sizeof(char) * 8)) - 1;
	int pixel_mask = (1 << (sizeof(char) * 16)) - 1;
	int col_pixels;
	int pixel_i, pixel_j;
	int palette_type;
	int dummy;
	int i, n;

	unsigned int pixel_id;
	unsigned int col_data[3];

	unsigned char pixel_r[4], pixel_g[4], pixel_b[4];

	char partial_image_name[256], input_image_name[256], output_image_name[256];

	char *pressure_unit = "(mmHg)";
	char *velocity_unit = "(m/s)";
	char *stress_unit = "(Pa)";

	glClear (GL_COLOR_BUFFER_BIT);

	glPointSize(1.F);
	glBegin (GL_POINTS);

	temp_file_ptr = fopen("temp_file.txt", "r");

	for (i = -1; i < image_count; i++) {
		fscanf(temp_file_ptr, "%s\n", partial_image_name);
	}
	fclose(temp_file_ptr);

	sprintf(input_image_name, "%s/%s", input_path, partial_image_name);

	std::string lImageNameNoExtension = std::string(partial_image_name).substr(
			0, std::string(partial_image_name).rfind('.'));

	sprintf(output_image_name, "%s/%s.ppm", output_path,
			lImageNameNoExtension.c_str());

	input_image_file = fopen(input_image_name, "r");
	xdrstdio_create(&xdr_image_file, input_image_file, XDR_DECODE);

	xdr_int(&xdr_image_file, &mode);
	xdr_float(&xdr_image_file, &pressure_threshold_min);
	xdr_float(&xdr_image_file, &pressure_threshold_max);
	xdr_float(&xdr_image_file, &velocity_threshold_max);
	xdr_float(&xdr_image_file, &stress_threshold_max);

	xdr_int(&xdr_image_file, &dummy);
	xdr_int(&xdr_image_file, &dummy);
	xdr_int(&xdr_image_file, &col_pixels);

	for (n = 0; n < col_pixels; n++) {
		xdr_u_int(&xdr_image_file, &pixel_id);

		xdr_u_int(&xdr_image_file, &col_data[0]);
		xdr_u_int(&xdr_image_file, &col_data[1]);
		xdr_u_int(&xdr_image_file, &col_data[2]);

		pixel_i = (pixel_id >> (2 * bits_per_char)) & pixel_mask;
		pixel_j = (pixel_id) & pixel_mask;

		pixel_r[0] = (col_data[0] >> (3 * bits_per_char)) & colour_mask;
		pixel_g[0] = (col_data[0] >> (2 * bits_per_char)) & colour_mask;
		pixel_b[0] = (col_data[0] >> (1 * bits_per_char)) & colour_mask;
		pixel_r[1] = (col_data[0] >> (0 * bits_per_char)) & colour_mask;
		pixel_g[1] = (col_data[1] >> (3 * bits_per_char)) & colour_mask;
		pixel_b[1] = (col_data[1] >> (2 * bits_per_char)) & colour_mask;
		pixel_r[2] = (col_data[1] >> (1 * bits_per_char)) & colour_mask;
		pixel_g[2] = (col_data[1] >> (0 * bits_per_char)) & colour_mask;
		pixel_b[2] = (col_data[2] >> (3 * bits_per_char)) & colour_mask;
		pixel_r[3] = (col_data[2] >> (2 * bits_per_char)) & colour_mask;
		pixel_g[3] = (col_data[2] >> (1 * bits_per_char)) & colour_mask;
		pixel_b[3] = (col_data[2] >> (0 * bits_per_char)) & colour_mask;

		glColor3f(pixel_r[0] * (1.F / 255.F), pixel_g[0] * (1.F / 255.F),
				pixel_b[0] * (1.F / 255.F));
		glVertex2f(pixel_i, pixel_j + (pixels_y >> 1));

		glColor3f(pixel_r[1] * (1.F / 255.F), pixel_g[1] * (1.F / 255.F),
				pixel_b[1] * (1.F / 255.F));
		glVertex2f(pixel_i + (pixels_x >> 1), pixel_j + (pixels_y >> 1));

		glColor3f(pixel_r[2] * (1.F / 255.F), pixel_g[2] * (1.F / 255.F),
				pixel_b[2] * (1.F / 255.F));
		glVertex2f(pixel_i, pixel_j);

		glColor3f(pixel_r[3] * (1.F / 255.F), pixel_g[3] * (1.F / 255.F),
				pixel_b[3] * (1.F / 255.F));

		glVertex2f(pixel_i + (pixels_x >> 1), pixel_j);
	}
	glEnd();

	xdr_destroy(&xdr_image_file);
	fclose(input_image_file);

	visVisualiseColourBar(RAINBOW, pixels_x * 0.1F, pixels_y * 0.5F + 20.F,
			pixels_x * 0.4F, pixels_y * 0.5F + 40.F);
	visVisualiseBarPars(velocity_unit, 0.0F, velocity_threshold_max,
			pixels_x * 0.1F, pixels_y * 0.5F + 20.F, pixels_x * 0.4F);

	visVisualiseColourBar(RAINBOW, pixels_x * 0.6F, pixels_y * 0.5F + 20.F,
			pixels_x * 0.9F, pixels_y * 0.5F + 40.F);
	visVisualiseBarPars(stress_unit, 0.0F, stress_threshold_max,
			pixels_x * 0.6F, pixels_y * 0.5F + 20.F, pixels_x * 0.9F);

	if (mode == 0) {
		palette_type = RAINBOW;
	} else if (mode == 1) {
		palette_type = HALF_RAINBOW;
	} else {
		palette_type = GREY;
	}
	visVisualiseColourBar(palette_type, pixels_x * 0.1F, 20.F, pixels_x * 0.4F,
			40.F);
	visVisualiseBarPars(pressure_unit, pressure_threshold_min,
			pressure_threshold_max, pixels_x * 0.1F, 20.F, pixels_x * 0.4F);

	visVisualiseColourBar(palette_type, pixels_x * 0.6F, 20.F, pixels_x * 0.9F,
			40.F);
	visVisualiseBarPars(stress_unit, 0.0F, stress_threshold_max,
			pixels_x * 0.6F, 20.F, pixels_x * 0.9F);

	glutSwapBuffers();

	if (cycle_id <= 2) {
		SaveWindowImage(pixels_x, pixels_y, output_image_name);
	}
	if (++image_count == images - 1) {
		++cycle_id;
		image_count = 0;
	}
}

void KeybordFunction(unsigned char key, int x, int y) {
	if (key == 'q') {
		free(image_buffer);
		free(image_data);

		exit(0);
	}
}

void Reshape(GLsizei w, GLsizei h) {
	//Projection (pixels_x, pixels_y);
}

void usage(char *progname) {
	printf(
			"Usage: %s <input directory path> <output directory path> <'x' if you want to view images> \n",
			progname);
}

int main(int argc, char *argv[]) {
	int required_args = 3;

	if (argc < required_args) {
		usage(argv[0]);
		exit(1);
	}

	input_path = argv[1];
	output_path = argv[2];

	FILE *image_file;
	FILE *temp_file_ptr;
	XDR xdr_image_file;

	int mode;

	float pressure_threshold_min, pressure_threshold_max;
	float velocity_threshold_max, stress_threshold_max;

	char partial_image_name[256], image_name[256];
	char *first_command_part, *last_command_part;
	char my_command[256];

	first_command_part = "ls -l";
	last_command_part = "| wc > temp_file.txt";

	sprintf(my_command, "%s %s %s", first_command_part, input_path,
			last_command_part);
	system(my_command);

	temp_file_ptr = fopen("temp_file.txt", "r");
	fscanf(temp_file_ptr, "%i ", &images);
	fclose(temp_file_ptr);

	temp_file_ptr = fopen("temp_file.txt", "r");

	first_command_part = "ls -1";
	last_command_part = "> temp_file.txt";

	sprintf(my_command, "%s %s %s", first_command_part, input_path,
			last_command_part);
	system(my_command);

	rewind(temp_file_ptr);

	fscanf(temp_file_ptr, "%s\n", partial_image_name);
	sprintf(image_name, "%s/%s", input_path, partial_image_name);

	fclose(temp_file_ptr);

	image_file = fopen(image_name, "r");
	xdrstdio_create(&xdr_image_file, image_file, XDR_DECODE);

	xdr_int(&xdr_image_file, &mode);
	xdr_float(&xdr_image_file, &pressure_threshold_min);
	xdr_float(&xdr_image_file, &pressure_threshold_max);
	xdr_float(&xdr_image_file, &velocity_threshold_max);
	xdr_float(&xdr_image_file, &stress_threshold_max);

	xdr_int(&xdr_image_file, &pixels_x);
	xdr_int(&xdr_image_file, &pixels_y);

	pixels_x *= 2;
	pixels_y *= 2;

	xdr_destroy(&xdr_image_file);
	fclose(image_file);

	if(argc == 4 && argv[3][0] == 'x' && argv[3][1] == '\0')
	{
		glutInit(&argc, argv);
		OpenWindow(pixels_x, pixels_y);

		Projection(pixels_x, pixels_y);

		image_data = (unsigned char *) malloc(
		sizeof(unsigned char) * pixels_x * pixels_y * 3);

		image_buffer = (unsigned char *) malloc(
		sizeof(unsigned char) * pixels_x * 3);

		//glutReshapeFunc (Reshape);
		glutIdleFunc(Display);
		glutDisplayFunc(Display);
		glutKeyboardFunc(KeybordFunction);
		glutMainLoop();
	}

	return (0);
}
