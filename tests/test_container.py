import sparpy

def test_container():
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

    p = sparpy.Particle2
    p.position = [5,6]
    particles[3] = p
    assert particles[3].position[0] == 5
    assert particles[3].position[1] == 6


