#include <cstdio>
#include <cstdlib>
#include <cmath>



bool write_matrix_to_file(const char * filename, const double * matrix, size_t num_rows, size_t num_cols)
{
    FILE * file = fopen(filename, "wb");
    if(file == nullptr)
    {
        fprintf(stderr, "Cannot open output file\n");
        return false;
    }

    fwrite(&num_rows, sizeof(size_t), 1, file);
    fwrite(&num_cols, sizeof(size_t), 1, file);
    fwrite(matrix, sizeof(double), num_rows * num_cols, file);

    fclose(file);

    return true;
}



void set_initial_solution(double * heat, size_t nx, size_t ny, double bc_north, double bc_south, double bc_west, double bc_east)
{
    // boundary conditions
    for(size_t x = 1; x < nx-1; x++) heat[(ny-1) * nx + x     ] = bc_north;
    for(size_t x = 1; x < nx-1; x++) heat[     0 * nx + x     ] = bc_south;
    for(size_t y = 1; y < ny-1; y++) heat[     y * nx + 0     ] = bc_west;
    for(size_t y = 1; y < ny-1; y++) heat[     y * nx + (nx-1)] = bc_east;
    heat[     0 * nx + 0     ] = (bc_south + bc_west) / 2;
    heat[(ny-1) * nx + 0     ] = (bc_north + bc_west) / 2;
    heat[     0 * nx + (ny-1)] = (bc_south + bc_east) / 2;
    heat[(ny-1) * nx + (ny-1)] = (bc_north + bc_east) / 2;

    // initial interior values - average of the boundary
    double initial_val = ((nx-1)*bc_north + (nx-1)*bc_south + (ny-1)*bc_west + (ny-1)*bc_east) / (2*nx + 2*ny - 4);
    for(size_t y = 1; y < ny-1; y++)
    {
        for(size_t x = 1; x < nx-1; x++)
        {
            heat[y * nx + x] = initial_val;
        }
    }
}



void copy_boundary(const double * source, double * destination, size_t nx, size_t ny)
{
    for(size_t x = 0; x < nx;   x++) destination[(ny-1) * nx + x     ] = source[(ny-1) * nx + x     ];
    for(size_t x = 0; x < nx;   x++) destination[     0 * nx + x     ] = source[     0 * nx + x     ];
    for(size_t y = 1; y < ny-1; y++) destination[     y * nx + 0     ] = source[     y * nx + 0     ];
    for(size_t y = 1; y < ny-1; y++) destination[     y * nx + (nx-1)] = source[     y * nx + (nx-1)];
}



void copy_heat(const double * source, double * destination, size_t nx, size_t ny)
{
    for(size_t y = 0; y < ny; y++)
    {
        for(size_t x = 0; x < nx; x++)
        {
            destination[y * nx + x] = source[y * nx + x];
        }
    }
}



void heat_iteration(const double * heat_curr, double * heat_next, size_t nx, size_t ny)
{
    for(size_t y = 1; y < ny-1; y++)
    {
        for(size_t x = 1; x < nx-1; x++)
        {
            double north_val = heat_curr[(y+1) * nx + x];
            double south_val = heat_curr[(y-1) * nx + x];
            double west_val =  heat_curr[    y * nx + (x-1)];
            double east_val =  heat_curr[    y * nx + (x+1)];
            double new_val = (north_val + south_val + west_val + east_val) / 4.0;
            heat_next[y * nx + x] = new_val;
        }
    }
}



double calculate_max_diff(const double * heat_curr, const double * heat_next, size_t nx, size_t ny)
{
    double max_diff = 0.0;

    for(size_t y = 1; y < ny-1; y++)
    {
        for(size_t x = 1; x < nx-1; x++)
        {
            size_t idx = y * nx + x;
            double diff = std::abs(heat_next[idx] - heat_curr[idx]);
            if(diff > max_diff)
            {
                max_diff = diff;
            }
        }
    }

    return max_diff;
}



void solve_heat(double * heat, size_t nx, size_t ny, int max_iterations, double epsilon)
{
    double * heat_help = new double[nx * ny];
    copy_boundary(heat, heat_help, nx, ny);

    double max_diff;

    int num_iters;
    for(num_iters = 1; num_iters <= max_iterations; num_iters++)
    {
        double * heat_curr = ((num_iters % 2 == 0) ? heat_help : heat);
        double * heat_next = ((num_iters % 2 == 0) ? heat : heat_help);

        heat_iteration(heat_curr, heat_next, nx, ny);
        max_diff = calculate_max_diff(heat_curr, heat_next, nx, ny);
        if(max_diff < epsilon) break;
    }

    if(num_iters % 2 != 0)
    {
        copy_heat(heat_help, heat, nx, ny);
    }

    if(num_iters <= max_iterations)
    {
        printf("Iterations converged in %d iterations with max_diff=%e\n", num_iters, max_diff);
    }
    else
    {
        printf("Iterations did not converge in %d iterations, max_diff=%e\n", max_iterations, max_diff);
    }

    delete[] heat_help;
}





int main(int argc, char ** argv)
{
    printf("Usage: ./random_matrix nx ny output_file.bin max_iters\n");
    printf("All parameters are optional and have default values\n");
    printf("\n");

    size_t nx = 10;
    size_t ny = 10;
    const char * output_file = "io/heat.bin";
    int max_iterations = 1000000;
    double epsilon = 1e-3;
    double bc_north = 0.0;
    double bc_south = 100.0;
    double bc_west = 100.0;
    double bc_east = 100.0;

    if(argc > 1) nx = static_cast<size_t>(atoll(argv[1]));
    if(argc > 2) ny = static_cast<size_t>(atoll(argv[2]));
    if(argc > 3) output_file = argv[3];
    if(argc > 4) max_iterations = atoi(argv[4]);

    printf("Command line arguments:\n");
    printf("  nx:             %zu\n", nx);
    printf("  ny:             %zu\n", ny);
    printf("  output_file:    %s\n", output_file);
    printf("  max_iterations: %d\n", max_iterations);
    printf("\n");

    if((ssize_t)nx <= 0 || (ssize_t)ny <= 0 || max_iterations < 0)
    {
        fprintf(stderr, "Wrong argument value\n");
        return 1;
    }



    printf("Initializing the rectangle ...\n");
    double * heat = new double[nx * ny];
    set_initial_solution(heat, nx, ny, bc_north, bc_south, bc_west, bc_east);
    printf("Done\n");
    printf("\n");


    printf("Solving the heat equation ...\n");
    solve_heat(heat, nx, ny, max_iterations, epsilon);
    printf("Done\n");
    printf("\n");

    printf("Writing matrix to file ...\n");
    bool success_write = write_matrix_to_file(output_file, heat, ny, nx);
    if(!success_write)
    {
        fprintf(stderr, "Failed to save matrix\n");
        return 2;
    }
    printf("Done\n");
    printf("\n");

    delete[] heat;

    printf("Finished successfully\n");

    return 0;
}
