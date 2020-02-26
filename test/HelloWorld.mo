// HelloWorld: solve a trivial differential equation dx/dt = -a*x
// DE solution: x(t) = K * exp{-a*x}
class HelloWorld
  Real x(start=1);
  // parameter: constant during the simulation, but can be easily changed
  // We could rerun the simulation for another a without compilation
  parameter Real a = 1;
  equation
    // The equality sign = specifies a equation, NOT AN ASSIGNMENT
    der(x) = -a*x;
end HelloWorld;