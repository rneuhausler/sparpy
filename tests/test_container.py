import sparpy

def test_container_2d():
    particles = sparpy.Particles2(10)
    particles[0].position = [1,2]
    print particles[0].position
    assert particles[0].position[0] == 1
    assert particles[0].position[1] == 2

    p = particles[1]
    p.position = [3,4]
    assert particles[1].position[0] == 3
    assert particles[1].position[1] == 4

    particles[2] = p
    assert particles[2].position[0] == 3
    assert particles[2].position[1] == 4

    p = sparpy.Particle2()
    p.position = [5,6]
    particles[3] = p
    assert particles[3].position[0] == 5
    assert particles[3].position[1] == 6

    particles2 = particles
    particles2[0].position = [6,7]
    assert particles[0].position[0] == 6
    assert particles[0].position[1] == 7
    assert particles2[0].position[0] == 6
    assert particles2[0].position[1] == 7

def test_container_3d():
    particles = sparpy.Particles3(10)
    particles[0].position = [1,2,3]
    assert particles[0].position[0] == 1
    assert particles[0].position[1] == 2
    assert particles[0].position[2] == 3

    p = sparpy.Particle3()
    p.position = [5,6,7]
    particles[0] = p
    assert particles[0].position[0] == 5
    assert particles[0].position[1] == 6
    assert particles[0].position[2] == 7

def test_container_1d():
    particles = sparpy.Particles1(10)
    particles[0].position = [1]
    assert particles[0].position[0] == 1

    p = sparpy.Particle1()
    p.position = [5]
    particles[0] = p
    assert particles[0].position[0] == 5

def test_append_container():
    particles = sparpy.Particles2()
    p = sparpy.Particle2()
    p.position = [1,2]
    particles.append(p)
    assert particles[0].position[0] == 1
    assert particles[0].position[1] == 2


def test_print_particle():
    p = sparpy.Particle2()
    p.position = [1,2]
    print 'test_print_particle: ',p


