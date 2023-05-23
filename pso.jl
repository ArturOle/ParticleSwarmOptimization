using Plots
# using BenchmarkTools
using Base.Threads
plotly()


mutable struct Particle
    coord::Vector{Float64}
    value_at_point::Float64

    hist_best_location::Vector{Float64}
    hist_best_value::Float64

    move_vector::Vector

    Particle(pos_vect; search="max") = new(pos_vect, -Inf, [0,0], -Inf, zeros(size(pos_vect)))
    Particle(pos_vect, search="min") = new(pos_vect, Inf, [0,0], Inf, zeros(size(pos_vect)))
end


function math_function(x)
    multiplication_result = 1
    sin_sum_result = 1
    for x_i in x
        multiplication_result *= x_i 
        sin_sum_result *= sin(x_i)
    end

    return sin_sum_result * sqrt(abs(multiplication_result))
end


function plot_function!(s, math_function)
    x, y = 0:0.1:10, 0:0.1:10
    z = Surface((x,y)->math_function([x,y]), x, y)
    surface!(s, x, y, z)
end


function plot_map!(p, math_function)
    x, y = 0:0.1:10, 0:0.1:10
    z = Surface((x,y)->math_function([x,y]), x, y)
    plot!(p, x, y, z, st = [:surface, :contourf])
end


function plot_particles(particle_matrix::Matrix)
    s = scatter(legend=false)
    
    scatter!(
        [particle.coord[1] for particle in particle_matrix], 
        [particle.coord[2] for particle in particle_matrix],
        color=colorant"black",
        xlabel="x",
        ylabel="y"
    )

    return s
end


function plot_particles!(plot, particle_matrix::Matrix)
    scatter!(
        plot,
        [particle.coord[1] for particle in particle_matrix], 
        [particle.coord[2] for particle in particle_matrix],
        [math_function(particle.coord) for particle in particle_matrix],
        color=colorant"black",
        markersize=2
    )
end


function plot_particles_hist!(plot, particle_matrix::Matrix)
    scatter!(
        plot,
        [particle.hist_best_location[1] for particle in particle_matrix], 
        [particle.hist_best_location[2] for particle in particle_matrix],
        [math_function(particle.hist_best_location) for particle in particle_matrix],
        color=colorant"black",
        markersize=2
    )
end


function plot_particles2d!(plot, particle_matrix::Matrix)
    scatter!(
        plot,
        [particle.coord[1] for particle in particle_matrix], 
        [particle.coord[2] for particle in particle_matrix],
        color=colorant"white",
        markersize=5
    )
end


function plot_particles_hist2d!(plot, particle_matrix::Matrix)
    scatter!(
        plot,
        [particle.hist_best_location[1] for particle in particle_matrix], 
        [particle.hist_best_location[2] for particle in particle_matrix],
        color=colorant"white",
        markersize=5
    )
end


function init_particle_matrix(N::Int, range::Tuple=((1.0,19.0), (1.0,19.0)); task="max")
    particle_matrix = Matrix{Particle}(undef, N, N)
    x_even_dist = [x for x in range[1][1]:(range[1][2]-range[1][1])/N:range[1][2]-(range[1][2]-range[1][1])/N]
    y_even_dist = [y for y in range[2][1]:(range[2][2]-range[2][1])/N:range[2][2]-(range[2][2]-range[2][1])/N]

    @threads for (index, vali) in collect(enumerate(x_even_dist))
        for (jndex, valj) in enumerate(y_even_dist)
            particle_matrix[jndex, index] = Particle([vali, valj]; search=task)
        end
    end

    return particle_matrix
end

# function init_particle_matrix(N::Int, range::Tuple=((1.0,10.0), (1.0,10.0)); dims::Int, task="max")
#     particle_matrix = Array{Particle}(undef, dims=[N for i in 1:dims])
#     even_dist_matrix = [
#         [x for x in range[i][1]:(range[i][2]-range[i][1])/N:range[i][2]-(range[i][2]-range[i][1])/N]
#         for i in 1:dims
#     ]
#     indexes_vector = Vector{Float64}(undef, dims)
#     positions_vector = Vector{Float64}(undef, dims)

#     @threads for dist_vect in even_dist_matrix
#         for (index, val) in enumerate(dist_vect)
            
#         end
#     end
# end

function init_particle_matrix4(N::Int, range::Tuple=((1.0,10.0), (1.0,10.0)); task="max")
    particle_matrix = Array{Particle}(undef, N, N, N, N)
    x_even_dist = [x for x in range[1][1]:(range[1][2]-range[1][1])/N:range[1][2]-(range[1][2]-range[1][1])/N]
    y_even_dist = [y for y in range[2][1]:(range[2][2]-range[2][1])/N:range[2][2]-(range[2][2]-range[2][1])/N]
    x_even_dist2 = [x for x in range[1][1]:(range[1][2]-range[1][1])/N:range[1][2]-(range[1][2]-range[1][1])/N]
    y_even_dist2 = [y for y in range[2][1]:(range[2][2]-range[2][1])/N:range[2][2]-(range[2][2]-range[2][1])/N]

    @threads for (index, vali) in collect(enumerate(x_even_dist))
        for (jndex, valj) in enumerate(y_even_dist)
            for (jndez, valz) in enumerate(x_even_dist2)
                for (jndea, vala) in enumerate(y_even_dist2)
                    particle_matrix[jndex, index, jndez, jndea] = Particle([vali, valj, valz, vala]; search=task)

                end
            end
        end
    end

    return particle_matrix
end


function eval_particles_max!(particle_matrix, math_function)
    glob_best_value = -Inf
    glob_best_pos = []
    for particle in particle_matrix
        particle.value_at_point = math_function(particle.coord)
        if isless(particle.hist_best_value, particle.value_at_point)
            particle.hist_best_location = particle.coord
            particle.hist_best_value = particle.value_at_point
        end
        if isless(glob_best_value, particle.value_at_point)
            glob_best_value = particle.value_at_point
            glob_best_pos = particle.coord
        end
    end
    return (glob_best_value, glob_best_pos)
end


function eval_particles!(particle_matrix, math_function; search=max)
    glob_best_value = -Inf
    glob_best_pos = []
    for particle in particle_matrix
        particle.value_at_point = math_function(particle.coord)
        if isless(particle.hist_best_value, particle.value_at_point)
            particle.hist_best_location = particle.coord
            particle.hist_best_value = particle.value_at_point
        end
        if isless(glob_best_value, particle.value_at_point)
            glob_best_value = particle.value_at_point
            glob_best_pos = particle.coord
        end
    end
    return (glob_best_value, glob_best_pos)
end


function eval_particles_min!(particle_matrix, math_function)
    glob_best_value = Inf
    glob_best_pos = []
    for particle in particle_matrix
        particle.value_at_point = math_function(particle.coord)
        if isless(particle.value_at_point, particle.hist_best_value)
            particle.hist_best_location = particle.coord
            particle.hist_best_value = particle.value_at_point
        end
        if isless(particle.value_at_point, glob_best_value)
            glob_best_value = particle.value_at_point
            glob_best_pos = particle.coord
        end
    end
    return (glob_best_value, glob_best_pos)
end


function eval_particles!(particle_matrix, math_function; search=min)
    glob_best_value = Inf
    glob_best_pos = []
    for particle in particle_matrix
        particle.value_at_point = math_function(particle.coord)
        if isless(particle.value_at_point, particle.hist_best_value)
            particle.hist_best_location = particle.coord
            particle.hist_best_value = particle.value_at_point
        end
        if isless(particle.value_at_point, glob_best_value)
            glob_best_value = particle.value_at_point
            glob_best_pos = particle.coord
        end
    end
    return (glob_best_value, glob_best_pos)
end


function eval_particle_min!(particle, math_function, glob_best_value, glob_best_pos)

    particle.value_at_point = math_function(particle.coord)
    if isless(particle.value_at_point, particle.hist_best_value)
        particle.hist_best_location = particle.coord
        particle.hist_best_value = particle.value_at_point
    end
    if isless(particle.value_at_point, glob_best_value)
        glob_best_value = particle.value_at_point
        glob_best_pos = particle.coord
    end
    return (glob_best_value, glob_best_pos)
end


function eval_particle_max!(particle, math_function, glob_best_value, glob_best_pos)

    particle.value_at_point = math_function(particle.coord)
    if particle.hist_best_value < particle.value_at_point
        particle.hist_best_location = particle.coord
        particle.hist_best_value = particle.value_at_point
    end
    if glob_best_value < particle.value_at_point
        glob_best_value = particle.value_at_point
        glob_best_pos = particle.coord
    end
    return (glob_best_value, glob_best_pos)
end


function calc_particle_next_v!(particle, glob_best_value, w_k, r_1, r_2, c_1=2, c_2=2)
    particle.move_vector = w_k*particle.move_vector .+
        c_1*r_1*(particle.hist_best_value - particle.value_at_point) .+
        c_2*r_2*(glob_best_value - particle.value_at_point)
end


function check_range_and_adjust(particle, range, index)
    if range[index][1] > particle.coord[index] 
        particle.coord[index] = range[index][1]
        particle.move_vector[index] = -particle.move_vector[index]
    elseif particle.coord[index] > range[index][2]
        particle.coord[index] = range[index][2]
        particle.move_vector[index] = -particle.move_vector[index]
    end
end


function calc_particle_next_x!(particle, range)
    particle.coord = particle.coord .+
    particle.move_vector./100

    for coord_index in eachindex(particle.coord)
        check_range_and_adjust(particle, range, coord_index)
    end
end


function pso(reps=1, math_function=math_function, range=((1.0,9.0),(1.0,9.0),(1.0,9.0),(1.0,9.0)), max_iter=300; searched="max")
    super_best_max = 0
    super_best_pos = []
    for i=1:reps
        if searched == "min"
            particle_matrix = init_particle_matrix(10, range; task="min")
            display(particle_matrix)
            (global_best_val, global_best_pos) = eval_particles_min!(particle_matrix, math_function)
            display(particle_matrix)
            display(global_best_pos)
            display(global_best_val)

            gr()
            s = plot(legend=false, xlabel="x", ylabel="y")
            plot_map!(s, math_function)
            plot_particles2d!(s, particle_matrix)
            display(s)

            plotly()
            s = scatter(legend=false, xlabel="x", ylabel="y")
            plot_function!(s, math_function)
            plot_particles!(s, particle_matrix)
            display(s)

            for iter in 1:max_iter
                (random_value_1, random_value_2) = rand(Float64, 2)
                wk = 0.9 - iter / (2*max_iter)
                for particle in particle_matrix
                    calc_particle_next_x!(particle, range)
                    calc_particle_next_v!(
                        particle,
                        global_best_val,
                        wk,
                        random_value_1,
                        random_value_2
                    )
                    
                    (global_best_val, global_best_pos) = eval_particle_min!(particle, math_function, global_best_val, global_best_pos)
                end
                s = plot(legend=false, xlabel="x", ylabel="y")
                plot_map!(s, math_function)
                plot_particles_hist2d!(s, particle_matrix)
                display(s)
                display(particle_matrix)

            end
            gr()
            s = plot(legend=false, xlabel="x", ylabel="y")
            plot_map!(s, math_function)
            plot_particles2d!(s, particle_matrix)


            plotly()
            s = scatter(legend=false, xlabel="x", ylabel="y")
            plot_function!(s, math_function)
            plot_particles!(s, particle_matrix)

            display((global_best_pos, global_best_val))
            if global_best_val > super_best_max
                super_best_max = global_best_val
                super_best_pos = global_best_pos
            end
        
        elseif searched == "max"
            particle_matrix = init_particle_matrix(5, range)
            (global_best_val, global_best_pos) = eval_particles_max!(particle_matrix, math_function)


            gr()
            s = plot(legend=false, xlabel="x", ylabel="y")
            plot_map!(s, math_function)
            plot_particles2d!(s, particle_matrix)
            display(s)

            # plotly()
            # s = scatter(legend=false, xlabel="x", ylabel="y")
            # plot_function!(s, math_function)
            # plot_particles!(s, particle_matrix)
            # display(s)


            anim = @animate for iter in 1:max_iter
                (random_value_1, random_value_2) = rand(Float64, 2)
                wk = 0.9 - iter / (2*max_iter)
                for particle in particle_matrix
                    calc_particle_next_x!(particle, range)
                    calc_particle_next_v!(
                        particle,
                        global_best_val,
                        wk,
                        random_value_1,
                        random_value_2
                    )
                    
                    (global_best_val, global_best_pos) = eval_particle_max!(particle, math_function, global_best_val, global_best_pos)
                end
                # gr()
                s = plot(legend=false, xlabel="x", ylabel="y")
                plot_map!(s, math_function)
                plot_particles_hist2d!(s, particle_matrix)
                # display(s)
                # display(particle_matrix)

            end

            display(gif(anim, fps=20))
            
            s = plot(legend=false, xlabel="x", ylabel="y")
            plot_map!(s, math_function)
            plot_particles_hist2d!(s, particle_matrix)
            display(s)
            display(particle_matrix)

            plotly()
            s = scatter(legend=false, xlabel="x", ylabel="y")
            plot_function!(s, math_function)
            plot_particles_hist!(s, particle_matrix)
            display(s)
            # display(particle_matrix)
            # display((global_best_pos, global_best_val))
            if global_best_val > super_best_max
                super_best_max = global_best_val
                super_best_pos = global_best_pos
            end
        end
    end
    display((super_best_pos, super_best_max))
end

pso()
