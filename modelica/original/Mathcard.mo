package Mathcard "Library of components to describe the Cardiovascular system and the cerebrospinal fluid"
  package Library "Library containing the blocks of the modelica mathcard lumped models library for cardiovascular and cerebrospinal system"
    package Vessels "Different kind of vessels"
      package I1O1
        package Autoregulating "This package contains windkessel element with Autoregulation possibilities. Each vessel is represented as a sequence of resistances, inertances and compliances through a linear system of differential equations."
          model CRL_LP_AVu "Implements the classical 0D model for vessels, obtained by linearizing the 1D model around the equilibrium state A=A0, P=0, Q=0. Vessel unstressed volume may vary through the addition of an external variation DeltaVu"
            extends Mathcard.Library.Connectors.I1O1;
            import Modelica.Constants.pi;
            Mathcard.Library.Connectors.RealInput Active_fes "fes-fesmin signal" annotation(
              Placement(visible = true, transformation(origin = {0.5088, 80.3926}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {-0.0, 40.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
            parameter Real R(unit = "mmHg.s/ml") = 1 "Compartment Resistance";
            parameter Real I(unit = "mmHg.s2/ml") = 1 "Compartment Inertance";
            parameter Real C(unit = "ml/mmHg") = 1 "Compartment Compliance";
            parameter Real V0 = 1 "Initial Compartment Volume";
            parameter Real Vu0 = 0 "Initial Compartment Unstressed Volume";
            parameter Real VuRef0 = 0 "Reference Compartment Unstressed Volume in absence of autoregulation";
            parameter Real AGain = 1 "Autoregulation signal gain";
            parameter Real ADelay = 1 "Autoregulation signal delay";
            parameter Real ATau = 1 "Autoregulation effect time constant";
          protected
            Real V(start = V0, unit = "ml") "Vessel Volume";
            Real PC(start = 1 / C * (V0 - Vu0), unit = "Pa") "Vessel Pressure";
            Real PCO(unit = "Pa");
            Real Qin(unit = "m3/s");
            Real Qout(unit = "m3/s");
            Real dV(unit = "m3/s");
            Real Vu(unit = "ml") "Compartment Resistance";
            Real sigmaVu(unit = "ml");
            Real deltaVu(start = Vu0 - VuRef0, unit = "ml");
          equation
            Inlet.P = PC;
            PCO = PC - Outlet.P;
            Qin = Inlet.Q;
            Qout = -Outlet.Q;
            dV = Qin - Qout;
            C * der(PC) = dV - der(Vu);
            der(V) = dV;
            if I == 0 then
              PCO = Qout * R;
            else
              I * der(Qout) = PCO - Qout * R;
            end if;
            sigmaVu = AGain * log(max(delay(Active_fes.Signal, ADelay), 0) + 1);
            der(deltaVu) = 1 / ATau * ((-deltaVu) + sigmaVu);
            Vu = VuRef0 + deltaVu;
            annotation(
              Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, origin = {30.0, 40.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "A", fontName = "Arial"), Rectangle(visible = true, fillColor = {204, 204, 204}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-70.0, -30.0}, {70.0, 30.0}})}));
          end CRL_LP_AVu;

          model CRL_LP_AR "Implements the classical 0D model for vessels, obtained by linearizing the 1D model around the equilibrium state A=A0, P=0, Q=0. Vessel resistance may vary through the addition of an external resistance variation DeltaR"
            extends Mathcard.Library.Connectors.I1O1;
            import Modelica.Constants.pi;
            Mathcard.Library.Connectors.RealInput Active_fes "fes-fesmin signal" annotation(
              Placement(visible = true, transformation(origin = {0.0, 81.4103}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {-0.0, 40.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
            parameter Real R0(unit = "mmHg.s2/ml") = 1 "Compartment Resistance at time 0";
            parameter Real RRef0(unit = "mmHg.s2/ml") = 1 "Compartment Resistance in absence of autoregulation stymulus";
            parameter Real I(unit = "mmHg.s2/ml") = 1 "Compartment Inertance";
            parameter Real C(unit = "ml/mmHg") = 1 "Compartment Compliance";
            parameter Real V0 = 1 "Initial Compartment Volume";
            parameter Real Vu0 = 0 "Initial Compartment Unstressed Volume";
            parameter Real AGain = 1 "Autoregulation signal gain";
            parameter Real ADelay = 1 "Autoregulation signal delay";
            parameter Real ATau = 1 "Autoregulation effect time constant";
          protected
            Real V(start = V0, unit = "ml") "Compartment Volume";
            Real PC(start = 1 / C * (V0 - Vu0), unit = "Pa") "Compartment Pressure";
            Real PCO(unit = "Pa");
            Real Qin(unit = "m3/s");
            Real Qout(unit = "m3/s");
            Real dV(unit = "m3/s");
            Real R(unit = "mmHg.s2/ml") "Compartment Resistance";
            Real sigmaR(unit = "mmHg.s2/ml");
            Real deltaR(start = R0 - RRef0, unit = "mmHg.s2/ml");
          equation
            Inlet.P = PC;
            PCO = PC - Outlet.P;
            Qin = Inlet.Q;
            Qout = -Outlet.Q;
            dV = Qin - Qout;
            C * der(PC) = dV;
            der(V) = dV;
            if I == 0 then
              PCO = Qout * R;
            else
              I * der(Qout) = PCO - Qout * R;
            end if;
            sigmaR = AGain * log(max(delay(Active_fes.Signal, ADelay), 0) + 1);
            der(deltaR) = 1 / ATau * ((-deltaR) + sigmaR);
            R = RRef0 + deltaR;
            annotation(
              Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, origin = {30.0, 40.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "A", fontName = "Arial"), Rectangle(visible = true, fillColor = {204, 204, 204}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-70.0, -30.0}, {70.0, 30.0}})}));
          end CRL_LP_AR;
          annotation(
            Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Autoregulation", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
            Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
            Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
        end Autoregulating;

        package Linear "This package contains classical windkessel element. Each vessel is represented as a sequence of resistances, inertances and compliances through a linear system of differential equations."
          model CRL_LP "Implements the classical 0D model for vessels, obtained by linearizing the 1D model around the equilibrium state A=A0, P=0, Q=0"
            extends Mathcard.Library.Connectors.I1O1;
            import Modelica.Constants.pi;
            parameter Real R(unit = "mmHg.s/ml") = 1 "Compartment Resistance";
            parameter Real I(unit = "mmHg.s2/ml") = 1 "Compartment Inertance";
            parameter Real C(unit = "ml/mmHg") = 1 "Compartment Compliance";
            parameter Real V0(unit = "ml") = 1 "Initial Compartment Volume";
            parameter Real Vu0 = 0 "Initial Compartment Unstressed Volume";
          protected
            Real V(start = V0, unit = "ml") "Vessel Volume";
            Real PC(start = 1 / C * (V0 - Vu0), unit = "Pa") "Vessel Pressure";
            Real PCO(unit = "Pa");
            Real Qin(unit = "m3/s");
            Real Qout(unit = "m3/s");
            Real dV(unit = "m3/s");
          equation
            Inlet.P = PC;
            PCO = PC - Outlet.P;
            Qin = Inlet.Q;
            Qout = -Outlet.Q;
            dV = Qin - Qout;
            C * der(PC) = dV;
            der(V) = dV;
            if I == 0 then
              PCO = Qout * R;
            else
              I * der(Qout) = PCO - Qout * R;
            end if;
            annotation(
              Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {204, 204, 204}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-70.0, -30.0}, {70.0, 30.0}})}));
          end CRL_LP;
          annotation(
            Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Linear", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
            Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
            Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
        end Linear;
        annotation(
          Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "1I 10", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
          Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
      end I1O1;

      package I1O1E1
        package Linear "This package contains classical windkessel element. Each vessel is represented as a sequence of resistances, inertances and compliances through a linear system of differential equations."
          model I1O1E1_nRLCLR_FP "Implements the classical 0D model for vessels, obtained by linearizing the 1D model around the equilibrium state A=A0, P=0, Q=0"
            extends Mathcard.Library.Connectors.I1O1E1;
            import Mathcard.Library.Vessels.Linear.*;
            parameter Real CL(unit = "cm") = 1 "Vessel Length";
            parameter Real Cru0(unit = "cm") = 1 "Vessel unstressed reference radius";
            parameter Real Chu0(unit = "cm") = 1 "Vessel unstressed reference thickness";
            parameter Real CEu0(unit = "MPa") = 1 "Vessel unstressed reference elastic modulus";
            parameter Real CV0(unit = "ml") = 1 "Initial Vessel Volume";
            parameter Real Cnu(unit = "1") = 1 "Vessel Poisson Ratio";
            parameter Real Crho(unit = "kg/m3") = 1060 "Blood Density";
            parameter Real Ceta(unit = "Pa.s") = 3.5 * 0.001 "Blood Dynamic Viscosity";
            parameter Integer N(final min = 1) = 1 "Number of lumped segments";
            Mathcard.Library.Vessels.I1O1E1.Linear.RLCLR_FP Segment[N + 1](each L = CL / (N + 1), each ru0 = Cru0, each Eu0 = CEu0, each nu = Cnu, each hu0 = Chu0, each rho = Crho, each eta = Ceta, each V0 = CV0 / (N + 1));
          equation
            connect(Inlet, Segment[1].Inlet);
            for i in 1:N loop
              connect(Segment[i].Outlet, Segment[i + 1].Inlet);
              connect(Segment[i].External, External);
            end for;
            connect(Segment[N + 1].Outlet, Outlet);
            connect(Segment[N + 1].External, External);
            annotation(
              Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {204, 204, 204}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-70.0, -30.0}, {70.0, 30.0}})}));
          end I1O1E1_nRLCLR_FP;

          model RLCLR_LP "Implements the classical 0D model for vessels, obtained by linearizing the 1D model around the equilibrium state A=A0, P=0, Q=0"
            extends Mathcard.Library.Connectors.I1O1E1;
            import Modelica.Constants.pi;
            parameter Real R(unit = "mmHg.s/ml") = 1 "Compartment Resistance";
            parameter Real I(unit = "mmHg.s2/ml") = 1 "Compartment Inertance";
            parameter Real C(unit = "ml/mmHg") = 1 "Compartment Compliance";
            parameter Real V0(unit = "ml") = 1 "Initial Compartment Volume";
            parameter Real Vu0 = 0 "Initial Compartment Unstressed Volume";
          protected
            Real V(start = V0, unit = "ml") "Vessel Volume";
            Real PTM(start = 1 / C * (V0 - Vu0), unit = "Pa") "Vessel Pressure";
            Real PC(unit = "Pa");
            Real PIC(unit = "Pa");
            Real PCO(unit = "Pa");
            Real Qin(unit = "m3/s");
            Real Qout(unit = "m3/s");
            Real dV(unit = "m3/s");
          equation
            PIC = Inlet.P - PC;
            PCO = PC - Outlet.P;
            PTM = PC - External.P;
            Qin = Inlet.Q;
            Qout = -Outlet.Q;
            dV = -External.dV;
            dV = Qin - Qout;
            C * der(PTM) = dV;
            der(V) = dV;
            if I == 0 then
              0 = PIC - Qin * R / 2;
              0 = PCO - Qout * R / 2;
            else
              I / 2 * der(Qin) = PIC - Qin * R / 2;
              I / 2 * der(Qout) = PCO - Qout * R / 2;
            end if;
            annotation(
              Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {204, 204, 204}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-70.0, -30.0}, {70.0, 30.0}})}));
          end RLCLR_LP;

          model RLCLR_FP "Implements the classical 0D model for vessels, obtained by linearizing the 1D model around the equilibrium state A=A0, P=0, Q=0"
            extends Mathcard.Library.Connectors.I1O1E1;
            import Modelica.Constants.pi;
            import Mathcard.Library.Quantities_and_units.*;
            parameter Real L(unit = "cm") = 1 "Vessel Length";
            parameter Real ru0(unit = "cm") = 1 "Vessel unstressed reference radius";
            parameter Real hu0(unit = "cm") = 1 "Vessel unstressed reference thickness";
            parameter Real Eu0(unit = "MPa") = 1 "Vessel unstressed reference elastic modulus";
            parameter Real V0(unit = "ml") = 1 "Initial Vessel Volume";
            parameter Real nu(unit = "1") = 1 "Vessel Poisson Ratio";
            parameter Real rho(unit = "kg/m3") = 1060 "Blood Density";
            parameter Real eta(unit = "Pa.s") = 3.5 * 0.001 "Blood Dynamic Viscosity";
          protected
            parameter Real R(unit = "mmHg.s/ml") = 8 * DynamicViscosityConversion(eta) * L / (Modelica.Constants.pi * ru0 ^ 4) "Vessel Resistance";
            parameter Real I(unit = "mmHg.s2/ml") = DensityConversion(rho) * L / (Modelica.Constants.pi * ru0 ^ 2) "Vessel Inertance";
            parameter Real C(unit = "ml/mmHg") = 2 * L * ru0 ^ 2 * (ru0 + hu0) ^ 2 * (1 + nu) * ElasticModulusConversion(Eu0) * hu0 * (hu0 + 2 * ru0) / (ElasticModulusConversion(Eu0) ^ 2 * hu0 ^ 2 * (hu0 + 2 * ru0) ^ 2) "Vessel Compliance";
            parameter Real Vu0 = L * ru0 ^ 2 * pi "Initial Compartment Unstressed Volume";
            Real V(start = V0, unit = "ml") "Vessel Volume";
            Real r(start = ru0, unit = "cm") "Vessel current radius";
            Real PTM(start = 1 / C * (V0 - Vu0), unit = "Pa") "Vessel Pressure";
            Real PC(unit = "Pa");
            Real PIC(unit = "Pa");
            Real PCO(unit = "Pa");
            Real Qin(unit = "m3/s");
            Real Qout(unit = "m3/s");
            Real dV(unit = "m3/s");
          equation
            PIC = Inlet.P - PC;
            PCO = PC - Outlet.P;
            PTM = PC - External.P;
            Qin = Inlet.Q;
            Qout = -Outlet.Q;
            dV = -External.dV;
            dV = Qin - Qout;
            C * der(PTM) = dV;
            der(V) = dV;
            if I == 0 then
              0 = PIC - Qin * R / 2;
              0 = PCO - Qout * R / 2;
            else
              I / 2 * der(Qin) = PIC - Qin * R / 2;
              I / 2 * der(Qout) = PCO - Qout * R / 2;
            end if;
            r = sqrt(V / (Modelica.Constants.pi * L));
            annotation(
              Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {204, 204, 204}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-70.0, -30.0}, {70.0, 30.0}})}));
          end RLCLR_FP;

          model nRLCLR_LP "Implements the classical 0D model for vessels, obtained by linearizing the 1D model around the equilibrium state A=A0, P=0, Q=0"
            extends Mathcard.Library.Connectors.I1O1E1;
            import Mathcard.Library.Vessels.Linear.*;
            parameter Real CR(unit = "mmHg.s/ml") = 1 "Compartment Resistance";
            parameter Real CI(unit = "mmHg.s2/ml") = 1 "Compartment Inertance";
            parameter Real CC(unit = "ml/mmHg") = 1 "Compartment Compliance";
            parameter Real CV0(unit = "ml") = 1 "Initial Compartment Volume";
            parameter Real CVu0 = 0 "Initial Compartment Unstressed Volume";
            parameter Integer N(final min = 1) = 1 "Number of lumped segments";
            Mathcard.Library.Vessels.I1O1E1.Linear.RLCLR_LP Segment[N + 1](each R = CR, each I = CI, each C = CC, each V0 = CV0, each Vu0 = CVu0);
          equation
            connect(Inlet, Segment[1].Inlet);
            for i in 1:N loop
              connect(Segment[i].Outlet, Segment[i + 1].Inlet);
              connect(Segment[i].External, External);
            end for;
            connect(Segment[N + 1].Outlet, Outlet);
            connect(Segment[N + 1].External, External);
            annotation(
              Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {204, 204, 204}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-70.0, -30.0}, {70.0, 30.0}})}));
          end nRLCLR_LP;
          annotation(
            Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Linear", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
            Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
            Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
        end Linear;
        annotation(
          Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "1I 10", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
          Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
      end I1O1E1;
      annotation(
        Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Vessels", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
        Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
        Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
        Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, origin = {-0.0, 45.0}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-85.0, -85.0}, {65.0, 35.0}}, textString = "Vessels", fontName = "Arial")}));
    end Vessels;

    package Quantities_and_units "Definition of quantities and units for the mathcard modelica library"
      type Pressure = Real(final quantity = "Pressure", final unit = "mmHg") annotation(
        Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
        Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
      type FlowRate = Real(final quantity = "FlowRate", final unit = "ml/s") annotation(
        Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
        Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
      type VolumeChangeRate = Real(final quantity = "FlowRate", final unit = "ml/s") annotation(
        Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
        Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));

      function DensityConversion "Convert between the usual Kg/m3 value for the density to the correct value to have pressures in mmHg and fluxes in ml/s"
        input Real RhoIn "Kg/m3 value";
        output Real RhoOut "mmHg.s2/cm2 value";
      algorithm
        RhoOut := RhoIn * 7.501 * 1e-07;
        annotation(
          Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
          Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, extent = {{-100, 20}, {-20, 100}}, textString = "Kg/m3", fontName = "Arial"), Text(visible = true, extent = {{20, -100}, {100, -20}}, textString = "mmHg.s2/cm2", fontName = "Arial")}),
          Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
      end DensityConversion;

      function ElasticModulusConversion "Convert between the usual MPa value for the elastic modulus to the correct value to have pressures in mmHg and fluxes in ml/s"
        input Real EIn "MPa value";
        output Real EOut "mmHg value";
      algorithm
        EOut := EIn * 7501;
        annotation(
          Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
          Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, extent = {{-100, 20}, {-20, 100}}, textString = "MPa", fontName = "Arial"), Text(visible = true, extent = {{20, -100}, {100, -20}}, textString = "mmHg", fontName = "Arial")}),
          Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
      end ElasticModulusConversion;

      function DynamicViscosityConversion "Convert between the usual MPa value for the elastic modulus to the correct value to have pressures in mmHg and fluxes in ml/s"
        input Real etaIn "Pa.s value";
        output Real etaOut "mmHg.s value";
      algorithm
        etaOut := etaIn * 7.501 * 0.001;
        annotation(
          Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
          Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, extent = {{-100, 20}, {-20, 100}}, textString = "Pa.s", fontName = "Arial"), Text(visible = true, extent = {{20, -100}, {100, -20}}, textString = "mmHg.s", fontName = "Arial")}),
          Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
      end DynamicViscosityConversion;
      annotation(
        Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Quantities & Units", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
        Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
        Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
    end Quantities_and_units;

    package Connectors "Connectors and partial models for Cardiovascular and CSF systems"
      connector Orifice "Orifice of a cardiovascular compartment (vessels, organs)"
        import Mathcard.Library.Quantities_and_units.*;
        Pressure P "Hydrostatic pressure at the node (orifice)";
        flow FlowRate Q "Flow rate flowing through the node (orifice)";
        annotation(
          Diagram(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
          Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Ellipse(visible = true, fillColor = {76, 76, 76}, fillPattern = FillPattern.Solid, extent = {{-90.0, -90.0}, {90.0, 90.0}})}));
      end Orifice;

      connector Wall "Connector to describe volume balance between different compartment"
        import Mathcard.Library.Quantities_and_units.*;
        Pressure P "External pressure at the compartment wall";
        flow VolumeChangeRate dV "VolumeChangeRate, positive if compartment is increasing its volume";
        annotation(
          Diagram(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
          Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {76, 76, 76}, fillPattern = FillPattern.Solid, extent = {{-90.0, -90.0}, {90.0, 90.0}})}));
      end Wall;

      partial model I1O1E1 "Partial model defining 1 Inlet Orifice, 1 Outlet Orifice and 1 External Wall"
        import Mathcard.Library.Connectors.*;
        Orifice Inlet "Inlet node" annotation(
          Placement(visible = true, transformation(origin = {-79.375, 0.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {-80.0, 0.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Orifice Outlet "Outlet node" annotation(
          Placement(visible = true, transformation(origin = {80.3926, -0.5088}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {80.0, -0.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Connectors.Wall External "Node accounting for compartment deformation (external pressure and volume change rate)" annotation(
          Placement(visible = true, transformation(origin = {-0.5088, -30.5288}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {0.0, -40.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        annotation(
          Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, origin = {-80.0, 30.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "I", fontName = "Arial"), Text(visible = true, origin = {80.0, 30.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "O", fontName = "Arial"), Text(visible = true, origin = {30.0, -40.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "E", fontName = "Arial")}));
      end I1O1E1;

      partial model I1O1 "Partial model defining 1 Inlet Orifice, 1 Outlet Orifice and 1 External Wall"
        import Mathcard.Library.Connectors.*;
        Mathcard.Library.Connectors.Orifice Inlet "Inlet node" annotation(
          Placement(visible = true, transformation(origin = {-79.375, 0.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {-80.0, 0.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Connectors.Orifice Outlet "Outlet node" annotation(
          Placement(visible = true, transformation(origin = {80.3926, -0.5088}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {80.0, 0.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        annotation(
          Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, origin = {-80.0, 30.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "I", fontName = "Arial"), Text(visible = true, origin = {80.0, 30.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "O", fontName = "Arial")}));
      end I1O1;

      connector RealInput
        input Real Signal "'input Real' as connector" annotation(
          defaultComponentName = "u",
          Documentation(info = "<html>
<p>
Connector with one input signal of type Real.
</p>
</html>"),
          Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid, points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}})}),
          Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid, points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}), Text(visible = true, fillColor = {0, 0, 127}, extent = {{-120, 60}, {100, 105}}, textString = "%name", fontName = "Arial")}));
        annotation(
          Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Ellipse(visible = true, fillColor = {0, 128, 0}, fillPattern = FillPattern.Solid, extent = {{-90.0, -90.0}, {90.0, 90.0}})}));
      end RealInput;

      connector RealOutput
        output Real Signal "'input Real' as connector" annotation(
          defaultComponentName = "u",
          Documentation(info = "<html>
<p>
Connector with one input signal of type Real.
</p>
</html>"),
          Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid, points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}})}),
          Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid, points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}), Text(visible = true, fillColor = {0, 0, 127}, extent = {{-120, 60}, {100, 105}}, textString = "%name", fontName = "Arial")}));
        annotation(
          Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Ellipse(visible = true, fillColor = {255, 128, 0}, fillPattern = FillPattern.Solid, extent = {{-90.0, -90.0}, {90.0, 90.0}})}));
      end RealOutput;
      annotation(
        Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Connectors", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
        Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
        Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
    end Connectors;

    package Sources_and_references "Pressure and flow sources and references elements for the cardiovascular system"
      package Sources
        model Bleeding
          Mathcard.Library.Connectors.Orifice BleedingInlet annotation(
            Placement(visible = true, transformation(origin = {-0.7556, 1.0685}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {0.1813, 56.2871}, extent = {{-40.1813, -40.1813}, {40.1813, 40.1813}}, rotation = 0)));
          parameter Real QBleeding(unit = "ml/s") = 0;
          parameter Real TstartBleeding(unit = "s") = 0;
          parameter Real TendBleeding(unit = "s") = 5;
        equation
          BleedingInlet.Q = if time > TstartBleeding and time < TendBleeding then QBleeding else 0;
          annotation(
            Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
            Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Ellipse(visible = true, origin = {-0.0, 0.0}, fillColor = {128, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}), Ellipse(visible = true, origin = {-0.0, -30.0}, fillColor = {128, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}), Ellipse(visible = true, origin = {-0.0, -61.0176}, fillColor = {128, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -11.0176}, {10.0, 11.0176}}), Ellipse(visible = true, origin = {-0.0, -85.0}, fillColor = {128, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-75.6823, -5.0}, {75.6823, 5.0}})}));
        end Bleeding;
        annotation(
          Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Sources", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
          Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
          Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
      end Sources;

      package References
        model PressureReference_dV "Gives a pressure reference value"
          Mathcard.Library.Connectors.Wall Node annotation(
            Placement(visible = true, transformation(origin = {-0.4992, 88.8601}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {-0.0, 90.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          parameter Real Pref(unit = "mmHg") = 0 "Value for the reference pressure";
        equation
          Node.P = Pref;
          annotation(
            Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Line(visible = true, points = {{-60.0, 50.0}, {60.0, 50.0}}), Line(visible = true, points = {{-40.0, 30.0}, {40.0, 30.0}}), Line(visible = true, points = {{-20.0, 10.0}, {20.0, 10.0}}), Line(visible = true, points = {{0.0, 90.0}, {0.0, 50.0}}), Text(visible = true, fillPattern = FillPattern.Solid, extent = {{-144.0, -60.0}, {138.0, 0.0}}, textString = "P=%Pref", fontName = "Arial")}));
        end PressureReference_dV;

        model PressureReference_Q "Gives a pressure reference value"
          Mathcard.Library.Connectors.Orifice Node annotation(
            Placement(visible = true, transformation(origin = {-0.4992, 88.8601}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {-0.0, 90.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          parameter Real Pref(unit = "mmHg") = 0 "Value for the reference pressure";
        equation
          Node.P = Pref;
          annotation(
            Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Line(visible = true, points = {{-60.0, 50.0}, {60.0, 50.0}}), Line(visible = true, points = {{-40.0, 30.0}, {40.0, 30.0}}), Line(visible = true, points = {{-20.0, 10.0}, {20.0, 10.0}}), Line(visible = true, points = {{0.0, 90.0}, {0.0, 50.0}}), Text(visible = true, fillPattern = FillPattern.Solid, extent = {{-144.0, -60.0}, {138.0, 0.0}}, textString = "P=%Pref", fontName = "Arial")}));
        end PressureReference_Q;
        annotation(
          Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "References", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
          Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
          Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
      end References;
      annotation(
        Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Sources and references", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
        Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
        Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
    end Sources_and_references;

    package Heart "Different kind of vessels"
      package Valves "This package contains windkessel element with Autoregulation possibilities. Each vessel is represented as a sequence of resistances, inertances and compliances through a linear system of differential equations."
        model Valve_AV "Implements the Left Atrium Model"
          import Mathcard.Library.Connectors.*;
          import Modelica.Constants.pi;
          Mathcard.Library.Connectors.Orifice Inlet annotation(
            Placement(visible = true, transformation(origin = {-0.5088, 50.3726}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {0.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          Mathcard.Library.Connectors.Orifice Outlet annotation(
            Placement(visible = true, transformation(origin = {0.0, -49.8638}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {0.0, -80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          parameter Real R(unit = "mmHg.s/ml");
        protected
          Real Q(unit = "ml/s") "Flux through the Valve";
        equation
          Q = Inlet.Q;
          Q = -Outlet.Q;
          Q = if noEvent(Inlet.P > Outlet.P) then (Inlet.P - Outlet.P) / R else 0;
          annotation(
            Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, origin = {30.0, 80.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "I", fontName = "Arial"), Text(visible = true, origin = {30.0, -80.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "O", fontName = "Arial"), Rectangle(visible = true, fillColor = {255, 255, 255}, extent = {{-40.0, -70.0}, {40.0, 70.0}}), Line(visible = true, origin = {22.5, 0.0}, points = {{17.5, 10.0}, {9.5553, 10.0}, {-2.5, -0.0}, {-12.5, -20.0}}, color = {153, 153, 153}, thickness = 5), Line(visible = true, points = {{0.0, 70.0}, {0.0, -70.0}}), Line(visible = true, origin = {-22.5, 0.0}, points = {{-17.5, 10.0}, {-9.5553, 10.0}, {2.5, -0.0}, {12.5, -20.0}}, color = {153, 153, 153}, thickness = 5), Line(visible = true, origin = {-2.5, 5.0}, points = {{-7.5, 5.0}, {12.5, 5.0}, {2.5, -15.0}, {-7.5, 5.0}}, thickness = 1)}));
        end Valve_AV;

        model Valve_VC "Implements the Left Atrium Model"
          import Mathcard.Library.Connectors.*;
          import Modelica.Constants.pi;
          Mathcard.Library.Connectors.Orifice Inlet annotation(
            Placement(visible = true, transformation(origin = {-0.5088, 50.3726}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {0.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          Mathcard.Library.Connectors.Orifice Outlet annotation(
            Placement(visible = true, transformation(origin = {0.0, -49.8638}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {-0.0, -80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          parameter Real kR(unit = "mmHg.s/ml") "Valve Resistance parameter";
        protected
          Real Q(unit = "ml/s") "Flux through the Valve";
        equation
          Q = Inlet.Q;
          Q = -Outlet.Q;
          Q = if Inlet.P > Outlet.P then (Inlet.P - Outlet.P) / (kR * Inlet.P) else 0;
          annotation(
            Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {255, 255, 255}, extent = {{-40.0, -70.0}, {40.0, 70.0}}), Line(visible = true, origin = {22.5, 0.0}, points = {{17.5, 10.0}, {9.5553, 10.0}, {-2.5, -0.0}, {-12.5, -20.0}}, color = {153, 153, 153}, thickness = 5), Line(visible = true, points = {{0.0, 70.0}, {0.0, -70.0}}), Line(visible = true, origin = {-2.5, 5.0}, points = {{-7.5, 5.0}, {12.5, 5.0}, {2.5, -15.0}, {-7.5, 5.0}}, thickness = 1), Line(visible = true, origin = {-22.5, 0.0}, points = {{-17.5, 10.0}, {-9.5553, 10.0}, {2.5, -0.0}, {12.5, -20.0}}, color = {153, 153, 153}, thickness = 5), Text(visible = true, origin = {30.0, 80.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "I", fontName = "Arial"), Text(visible = true, origin = {30.0, -80.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "O", fontName = "Arial")}));
        end Valve_VC;
        annotation(
          Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Autoregulation", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
          Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
          Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
      end Valves;

      package AutoregulatingChambers "This package contains windkessel element with Autoregulation possibilities. Each vessel is represented as a sequence of resistances, inertances and compliances through a linear system of differential equations."
        model Atrium "Model for an Atrial Cardiac Chamber"
          import Mathcard.Library.Connectors.*;
          import Modelica.Constants.pi;
          Mathcard.Library.Connectors.Orifice Inlet annotation(
            Placement(visible = true, transformation(origin = {-0.5088, 50.3726}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {-0.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          Mathcard.Library.Connectors.Orifice Outlet annotation(
            Placement(visible = true, transformation(origin = {0.0, -49.8638}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {0.0, -80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          parameter Real C(unit = "ml/mmHg") "Atrial Compliance";
          parameter Real V0(unit = "ml") "Initial Atrial Volume";
          parameter Real Vu0 = 0 "Initial Atrial Unstressed Volume";
        protected
          Real V(start = V0, unit = "ml") "Atrial Volume";
          Real PC(start = 1 / C * (V0 - Vu0), unit = "mmHg") "Atrial Pressure";
          Real Qin(unit = "ml/s") "Inlet flux into the Atrium";
          Real Qout(unit = "ml/s") "Outlet flux from the Atrium";
        equation
          Qin = Inlet.Q;
          Inlet.P = PC;
          C * der(PC) = Qin - Qout;
          Qout = -Outlet.Q;
          Outlet.P = PC;
          der(V) = Qin - Qout;
          annotation(
            Diagram(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
            Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Ellipse(visible = true, origin = {0.0, 0.0}, lineColor = {128, 0, 0}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Sphere, extent = {{-70.0, -70.0}, {70.0, 70.0}}), Text(visible = true, origin = {30.0, 80.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "I", fontName = "Arial"), Text(visible = true, origin = {30.0, -80.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "O", fontName = "Arial")}));
        end Atrium;

        model Ventricle "Model for a Ventricular Cardiac Chamber"
          import Mathcard.Library.Connectors.*;
          import Modelica.Constants.pi;
          Mathcard.Library.Connectors.Orifice Inlet annotation(
            Placement(visible = true, transformation(origin = {-0.5088, 50.3726}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {-0.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          Mathcard.Library.Connectors.Orifice Outlet annotation(
            Placement(visible = true, transformation(origin = {0.0, -49.8638}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {0.0, -80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          Mathcard.Library.Connectors.RealInput Active_fes annotation(
            Placement(visible = true, transformation(origin = {30.5288, -19.8438}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {90.0, -30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          Mathcard.Library.Connectors.RealInput Active_fev annotation(
            Placement(visible = true, transformation(origin = {26.4583, 18.8261}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {90.0, 30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          parameter Real V0(unit = "ml") = 1 "Initial Ventricular Volume";
          parameter Real kR(unit = "s/ml") = 1 "Aortic Valve Parameter";
          parameter Real kE(unit = "1/ml") = 1 "Monoexponential relation parameter";
          parameter Real P0(unit = "mmHg") = 1 "Monoexponential relation parameter";
          parameter Real kSys(unit = "s2") = 1 "Systole parameter";
          parameter Real TSys0(unit = "s") = 1 "Systole reference duration";
          parameter Real AGain_Emax = 1;
          parameter Real ADelay_Emax = 1;
          parameter Real ATau_Emax = 1;
          parameter Real Emax0 = 1;
          parameter Real EmaxRef0 = 1;
          parameter Real ATau_Ts = 1;
          parameter Real ATau_Tv = 1;
          parameter Real AGain_Ts = -1;
          parameter Real AGain_Tv = 1;
          parameter Real ADelay_Ts = 1;
          parameter Real ADelay_Tv = 1;
          parameter Real TRef0 = 1;
        protected
          Real V(start = V0, unit = "ml") "Ventricular Volume";
          Real PV(unit = "mmHg") "Ventricular Chamber Pressure";
          Real PVmax(unit = "mmHg") "Ventricular Pressure without the dynamic part in proximity of the outlet valve";
          Real Qin(unit = "ml/s") "Inlet flux into the Ventriculum";
          Real Qout(unit = "ml/s") "Outlet flux from the Ventriculum";
          Real TSys(unit = "s") "Systole duration";
          Real xi(start = 0, unit = "1") "Auxiliary Variable";
          Real u(unit = "1") "Fraction of the cardiac cycle";
          Real phi(unit = "1") "Ventricle activation function";
          Real Emax;
          Real sigmaEmax;
          Real deltaEmax(start = Emax0 - EmaxRef0);
          Real sigmaTs;
          Real sigmaTv;
          Real deltaTs(start = 0);
          Real deltaTv(start = 0);
          Real HR;
        public 
          Real T;           // should be protected
        equation
          Qin = Inlet.Q;
          Qout = -Outlet.Q;
          PV = Inlet.P;
          PVmax = Outlet.P;
          der(V) = Qin - Qout;
          PV = PVmax * (1 - kR * Qout);
          PVmax = phi * Emax * (V - V0) + (1 - phi) * P0 * (exp(kE * V) - 1);
          Emax = deltaEmax + EmaxRef0;
          der(deltaEmax) = 1 / ATau_Emax * ((-deltaEmax) + sigmaEmax);
          sigmaEmax = AGain_Emax * log(max(delay(Active_fes.Signal, ADelay_Emax), 0) + 1);
          TSys = TSys0 - kSys / T;
          der(xi) = 1 / T;
          u = xi - floor(xi);
          phi = if 0 <= u and u * T <= TSys then sin(pi * T * u / TSys) ^ 2 else 0;
          sigmaTs = AGain_Ts * log(max(delay(Active_fes.Signal, ADelay_Ts), 0) + 1);
          der(deltaTs) = 1 / ATau_Ts * ((-deltaTs) + sigmaTs);
          sigmaTv = AGain_Tv * delay(Active_fev.Signal, ADelay_Tv);
          der(deltaTv) = 1 / ATau_Tv * ((-deltaTv) + sigmaTv);
          T = TRef0 + deltaTs + deltaTv;
          HR = 60 / T;
          annotation(
            Diagram(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
            Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, origin = {30.0, 80.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "I", fontName = "Arial"), Text(visible = true, origin = {30.0, -80.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "O", fontName = "Arial"), Text(visible = true, origin = {90.0, 50.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "V", fontName = "Arial"), Text(visible = true, origin = {90.0, -50.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "S", fontName = "Arial"), Rectangle(visible = true, lineColor = {128, 0, 0}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Sphere, extent = {{-80.0, -70.0}, {80.0, 70.0}})}));
        end Ventricle;
        annotation(
          Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Autoregulation", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
          Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
          Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
      end AutoregulatingChambers;
      annotation(
        Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Vessels", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
        Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
        Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
        Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, origin = {-0.0, 45.0}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-85.0, -85.0}, {65.0, 35.0}}, textString = "Vessels", fontName = "Arial")}));
    end Heart;

    package AutoregulationCenters "Different kind of vessels"
      model Autoregulation_Center
        import Mathcard.Library.Connectors.*;
        Mathcard.Library.Connectors.Orifice PressureProbe annotation(
          Placement(visible = true, transformation(origin = {-117.8719, 0.5343}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {-80.0, 0.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Connectors.RealOutput Active_fes annotation(
          Placement(visible = true, transformation(origin = {118.6275, 33.6581}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {80.0, 30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Connectors.RealOutput Active_fev annotation(
          Placement(visible = true, transformation(origin = {118.6275, -29.9183}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {80.0, -30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        parameter Real TauP = 2.076;
        parameter Real TauZ = 6.37;
        parameter Real Pn = 92;
        parameter Real ka = 11.758;
        parameter Real fmin = 2.52;
        parameter Real fmax = 47.78;
        parameter Real fesinf = 2.1;
        parameter Real fes0 = 16.11;
        parameter Real fesmin = 2.66;
        parameter Real kes = 0.0675;
        parameter Real fev0 = 3.2;
        parameter Real fevinf = 6.3;
        parameter Real fcs0 = 25;
        parameter Real kev = 7.06;
      protected
        Real Qin;
        Real Pcs "Pressure at the carotid sinus";
        Real Ptilde "Auxiliary autoregulating pressure";
        Real fcs;
        Real fes;
      equation
        PressureProbe.Q = Qin;
        PressureProbe.P = Pcs;
        Qin = 0;
        TauP * der(Ptilde) = Pcs + TauZ * der(Pcs) - Ptilde;
        fcs = (fmin + fmax * exp((Ptilde - Pn) / ka)) / (1 + exp((Ptilde - Pn) / ka));
        fes = fesinf + (fes0 - fesinf) * exp(-kes * fcs);
        Active_fev.Signal = (fev0 + fevinf * exp((fcs - fcs0) / kev)) / (1 + exp((fcs - fcs0) / kev));
        Active_fes.Signal = fes - fesmin;
//  annotation(Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})), Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {64, 128, 0}, fillPattern = FillPattern.Solid, extent = {{-70.0, -40.0}, {70.0, 40.0}}), Text(visible = true, origin = {-60.0, -0.0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "Pc", fontSize = 50, fontName = "Arial"), Text(visible = true, origin = {55.0, 30.0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-15.0, -10.0}, {15.0, 10.0}}, textString = "Afes", fontSize = 50, fontName = "Arial"), Text(visible = true, origin = {55.0, -30.0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-15.0, -10.0}, {15.0, 10.0}}, textString = "Afev", fontSize = 50, fontName = "Arial"), Text(visible = true, origin = {-0.0, -50.0}, fillColor = {0, 128, 0}, fillPattern = FillPattern.Solid, extent = {{-70.0, -10.0}, {70.0, 10.0}}, textString = "Autoregulation Center", fontSize = 50, fontName = "Arial")}));
      end Autoregulation_Center;
      annotation(
        Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Vessels", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
        Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
        Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
        Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, origin = {-0.0, 45.0}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-85.0, -85.0}, {65.0, 35.0}}, textString = "Vessels", fontName = "Arial")}));
    end AutoregulationCenters;

    package Devices
      model VAD
        import Mathcard.Library.Connectors.*;
        parameter Real RPM;
        Mathcard.Library.Connectors.Orifice Inlet annotation(
          Placement(visible = true, transformation(origin = {0.0, 60.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {0.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Connectors.Orifice Outlet annotation(
          Placement(visible = true, transformation(origin = {0.0, -51.6731}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {0.0, -80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
      public
        // this should be protected
        Real Q;
        Real dP;
      equation
        Inlet.Q = Q;
        Outlet.Q = -Q;
        dP = Outlet.P - Inlet.P;
        Q = if RPM == 10000 then min(306.65 - 4.03 * dP + 0.0127 * dP ^ 2, 133) else if RPM == 9000 then min(245.16 - 3.6 * dP + 0.0119 * dP ^ 2, 116.667) else if RPM == 8000 then min(195.85 - 3.506 * dP + 0.0136 * dP ^ 2, 100) else 0;
        annotation(
          Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
          Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {255, 128, 0}, fillPattern = FillPattern.Solid, extent = {{-20.0, -70.0}, {20.0, 70.0}}), Text(visible = true, fillColor = {128, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-16.8042, -10.5833}, {16.8042, 10.5833}}, textString = "VAD", fontName = "Arial"), Text(visible = true, origin = {30.0, 80.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "I", fontName = "Arial"), Text(visible = true, origin = {30.0, -80.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "O", fontName = "Arial")}));
      end VAD;

      model ECMO
        import Mathcard.Library.Connectors.*;
        parameter Real RPM;
        Mathcard.Library.Connectors.Orifice Inlet annotation(
          Placement(visible = true, transformation(origin = {0.0, 60.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {0.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Connectors.Orifice Outlet annotation(
          Placement(visible = true, transformation(origin = {0.0, -51.6731}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {0.0, -80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
      protected
        Real Q;
        Real dP;
      equation
        Inlet.Q = Q;
        Outlet.Q = -Q;
        dP = Outlet.P - Inlet.P;
        Q = if RPM == 10000 then min(306.65 - 4.03 * dP + 0.0127 * dP ^ 2, 133) else if RPM == 9000 then min(245.16 - 3.6 * dP + 0.0119 * dP ^ 2, 116.667) else if RPM == 8000 then min(195.85 - 3.506 * dP + 0.0136 * dP ^ 2, 100) else 0;
        annotation(
          Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
          Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {255, 128, 0}, fillPattern = FillPattern.Solid, extent = {{-20.0, -70.0}, {20.0, 70.0}}), Text(visible = true, fillColor = {128, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-16.8042, -10.5833}, {16.8042, 10.5833}}, textString = "VAD", fontName = "Arial"), Text(visible = true, origin = {30.0, 80.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "I", fontName = "Arial"), Text(visible = true, origin = {30.0, -80.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "O", fontName = "Arial")}));
      end ECMO;

      model VAD2 "Reimplementation of VAD model for Heart Mate III, i.e. with a new HQ curve"
        import Mathcard.Library.Connectors.*;
        import Modelica.Blocks.Sources.Pulse;
        // Parameters for Artificial Pulse
        parameter Real RPM_c;
        parameter Real HMIII_Pulse_Period = 2.0;
        parameter Real HMIII_Pulse_Amplitude = 2000;
        parameter Real HMIII_Pulse_startTime = 0.0;
        parameter Real HMIII_Pulse_Up_Width = 10; // percentage of period
        parameter Real HMIII_Pulse_Down_Width = 7.5; // percentage of period
        // Pulse generators
        Modelica.Blocks.Sources.Pulse PulseGeneratorUp(period = HMIII_Pulse_Period, amplitude=HMIII_Pulse_Amplitude, offset=RPM_c, width=HMIII_Pulse_Up_Width, startTime = HMIII_Pulse_startTime + HMIII_Pulse_Period*HMIII_Pulse_Down_Width/100);
        Modelica.Blocks.Sources.Pulse PulseGeneratorDown(period= HMIII_Pulse_Period, amplitude=HMIII_Pulse_Amplitude, offset=0, width=HMIII_Pulse_Down_Width, startTime=HMIII_Pulse_Period);
        // Inlet/Outlet
        Mathcard.Library.Connectors.Orifice Inlet annotation(
          Placement(visible = true, transformation(origin = {0.0, 60.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {0.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Connectors.Orifice Outlet annotation(
          Placement(visible = true, transformation(origin = {0.0, -51.6731}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0), iconTransformation(origin = {0.0, -80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
      //public
        // this should be protected
      protected
        Real Q;
        Real dP;
      public
        Real RPM;
      equation
        // Use pulse function to implement artificial pulse
        RPM = PulseGeneratorUp.y - PulseGeneratorDown.y;
      //
        Inlet.Q = Q;
        Outlet.Q = -Q;
        dP = Outlet.P - Inlet.P;
        // For now we artificially prevent any backlow
        // by setting the min of Q to 0
        Q = max(0, if RPM == 3000 then if dP < 12.33 then 48.43 else (-0.12205 * dP ^ 2) + 3.0093 * dP + 29.886 else if RPM == 4000 then if dP < 18.94 then 75.27 else (-0.048803 * dP ^ 2) + 1.8484 * dP + 57.770 else if RPM == 5000 then if dP < 26.36 then 101.17 else (-0.024594 * dP ^ 2) + 1.29678 * dP + 84.081 else if RPM == 6000 then if dP < 33.71 then 120.83 else (-0.012153 * dP ^ 2) + 0.81934 * dP + 107.02 else if RPM == 7000 then if dP < 78.76 then 136.57 else (-0.013720 * dP ^ 2) + 2.1611 * dP + 51.474 else if RPM == 8000 then if dP < 63.73 then 166.28 else (-0.0055313 * dP ^ 2) + 0.70504 * dP + 143.71 else if dP < 93.82 then 178.95 else (-0.0042916 * dP ^ 2) + 0.80525 * dP + 141.18);
        //Q = 0;
        annotation(
          Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
          Icon(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {255, 128, 0}, fillPattern = FillPattern.Solid, extent = {{-20.0, -70.0}, {20.0, 70.0}}), Text(visible = true, fillColor = {128, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-16.8042, -10.5833}, {16.8042, 10.5833}}, textString = "VAD", fontName = "Arial"), Text(visible = true, origin = {30.0, 80.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "I", fontName = "Arial"), Text(visible = true, origin = {30.0, -80.0}, fillPattern = FillPattern.Solid, extent = {{-10.0, -10.0}, {10.0, 10.0}}, textString = "O", fontName = "Arial")}));
      end VAD2;
      annotation(
        Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Devices", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
        Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
    end Devices;
    annotation(
      Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Library", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
      Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
  end Library;

  package Applications "Package containing the application developed with the lumped models mathcard modelica library"
    package Ursino1998
      model Ursino1998Model
        import Mathcard.Library.*;
        extends Mathcard.Applications.Ursino1998.ModelParametersNH;
        Mathcard.Library.Heart.AutoregulatingChambers.Atrium LeftAtrium(C = Param_LeftAtrium_C, V0 = Param_LeftAtrium_V0, Vu0 = Param_LeftAtrium_Vu0) annotation(
          Placement(visible = true, transformation(origin = {30.0, 30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Heart.Valves.Valve_AV MitralicValve(R = Param_MitralicValve_R) annotation(
          Placement(visible = true, transformation(origin = {30.0, 10.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Heart.AutoregulatingChambers.Ventricle LeftVentricle(V0 = Param_LeftVentricle_V0, kR = Param_LeftVentricle_kR, Emax0 = Param_LeftVentricle_Emax0, EmaxRef0 = Param_LeftVentricle_EmaxRef0, AGain_Emax = Param_LeftVentricle_AGain_Emax, ADelay_Emax = Param_LeftVentricle_ADelay_Emax, ATau_Emax = Param_LeftVentricle_ATau_Emax, P0 = Param_LeftVentricle_P0, kE = Param_LeftVentricle_kE, TSys0 = Param_LeftVentricle_TSys0, kSys = Param_LeftVentricle_kSys, TRef0 = Param_LeftVentricle_TRef0, AGain_Ts = Param_LeftVentricle_AGain_Ts, ADelay_Ts = Param_LeftVentricle_ADelay_Ts, ATau_Ts = Param_LeftVentricle_ATau_Ts, AGain_Tv = Param_LeftVentricle_AGain_Tv, ADelay_Tv = Param_LeftVentricle_ADelay_Tv, ATau_Tv = Param_LeftVentricle_ATau_Tv) annotation(
          Placement(visible = true, transformation(origin = {30.0, -10.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Heart.Valves.Valve_VC AorticValve(kR = Param_AorticValve_kR) annotation(
          Placement(visible = true, transformation(origin = {30.0, -30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Vessels.I1O1.Linear.CRL_LP SystemicArteries(R = Param_SystemicArteries_R, I = Param_SystemicArteries_I, C = Param_SystemicArteries_C, V0 = Param_SystemicArteries_V0, Vu0 = Param_SystemicArteries_Vu0) annotation(
          Placement(visible = true, transformation(origin = {70.0, -80.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Vessels.I1O1.Autoregulating.CRL_LP_AR SplanchnicPeripheralCirculation(R0 = Param_SplanchnicPeripheralCirculation_R0, RRef0 = Param_SplanchnicPeripheralCirculation_RRef0, AGain = Param_SplanchnicPeripheralCirculation_AGain, ADelay = Param_SplanchnicPeripheralCirculation_ADelay, ATau = Param_SplanchnicPeripheralCirculation_ATau, I = Param_SplanchnicPeripheralCirculation_I, C = Param_SplanchnicPeripheralCirculation_C, V0 = Param_SplanchnicPeripheralCirculation_V0, Vu0 = Param_SplanchnicPeripheralCirculation_Vu0) annotation(
          Placement(visible = true, transformation(origin = {0.0, -90.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Vessels.I1O1.Autoregulating.CRL_LP_AR ExtraSplanchnicPeripheralCirculation(R0 = Param_ExtraSplanchnicPeripheralCirculation_R0, RRef0 = Param_ExtraSplanchnicPeripheralCirculation_RRef0, AGain = Param_ExtraSplanchnicPeripheralCirculation_AGain, ADelay = Param_ExtraSplanchnicPeripheralCirculation_ADelay, ATau = Param_ExtraSplanchnicPeripheralCirculation_ATau, I = Param_ExtraSplanchnicPeripheralCirculation_I, C = Param_ExtraSplanchnicPeripheralCirculation_C, V0 = Param_ExtraSplanchnicPeripheralCirculation_V0, Vu0 = Param_ExtraSplanchnicPeripheralCirculation_Vu0) annotation(
          Placement(visible = true, transformation(origin = {0.0, -70.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -180)));
        Mathcard.Library.Vessels.I1O1.Autoregulating.CRL_LP_AVu SplanchnicVeins(R = Param_SplanchnicVeins_R, I = Param_SplanchnicVeins_I, C = Param_SplanchnicVeins_C, V0 = Param_SplanchnicVeins_V0, Vu0 = Param_SplanchnicVeins_Vu0, VuRef0 = Param_SplanchnicVeins_VuRef0, AGain = Param_SplanchnicVeins_AGain, ADelay = Param_SplanchnicVeins_ADelay, ATau = Param_SplanchnicVeins_ATau) annotation(
          Placement(visible = true, transformation(origin = {-70.0, -90.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Vessels.I1O1.Autoregulating.CRL_LP_AVu ExtraSplanchnicVeins(R = Param_ExtraSplanchnicVeins_R, I = Param_ExtraSplanchnicVeins_I, C = Param_ExtraSplanchnicVeins_C, V0 = Param_ExtraSplanchnicVeins_V0, Vu0 = Param_ExtraSplanchnicVeins_Vu0, VuRef0 = Param_ExtraSplanchnicVeins_VuRef0, AGain = Param_ExtraSplanchnicVeins_AGain, ADelay = Param_ExtraSplanchnicVeins_ADelay, ATau = Param_ExtraSplanchnicVeins_ATau) annotation(
          Placement(visible = true, transformation(origin = {-70.0, -70.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -180)));
        Mathcard.Library.Heart.AutoregulatingChambers.Atrium RightAtrium(C = Param_RightAtrium_C, V0 = Param_RightAtrium_V0, Vu0 = Param_RightAtrium_Vu0) annotation(
          Placement(visible = true, transformation(origin = {-30.0, 30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Heart.Valves.Valve_AV TricuspidValve(R = Param_TricuspidValve_R) annotation(
          Placement(visible = true, transformation(origin = {-30.0, 10.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Heart.AutoregulatingChambers.Ventricle RightVentricle(V0 = Param_RightVentricle_V0, kR = Param_RightVentricle_kR, Emax0 = Param_RightVentricle_Emax0, EmaxRef0 = Param_RightVentricle_EmaxRef0, AGain_Emax = Param_RightVentricle_AGain_Emax, ADelay_Emax = Param_RightVentricle_ADelay_Emax, ATau_Emax = Param_RightVentricle_ATau_Emax, P0 = Param_RightVentricle_P0, kE = Param_RightVentricle_kE, TSys0 = Param_RightVentricle_TSys0, kSys = Param_RightVentricle_kSys, TRef0 = Param_RightVentricle_TRef0, AGain_Ts = Param_RightVentricle_AGain_Ts, ADelay_Ts = Param_RightVentricle_ADelay_Ts, ATau_Ts = Param_RightVentricle_ATau_Ts, AGain_Tv = Param_RightVentricle_AGain_Tv, ADelay_Tv = Param_RightVentricle_ADelay_Tv, ATau_Tv = Param_RightVentricle_ATau_Tv) annotation(
          Placement(visible = true, transformation(origin = {-30.0, -10.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Heart.Valves.Valve_VC PulmonaryValve(kR = Param_PulmonaryValve_kR) annotation(
          Placement(visible = true, transformation(origin = {-30.0, -30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Vessels.I1O1.Linear.CRL_LP PulmonaryArteries(R = Param_PulmonaryArteries_R, I = Param_PulmonaryArteries_I, C = Param_PulmonaryArteries_C, V0 = Param_PulmonaryArteries_V0, Vu0 = Param_PulmonaryArteries_Vu0) annotation(
          Placement(visible = true, transformation(origin = {-50.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -360)));
        Mathcard.Library.Vessels.I1O1.Linear.CRL_LP PulmonaryPeriphericalCirculation(R = Param_PulmonaryPeriphericalCirculation_R, I = Param_PulmonaryPeriphericalCirculation_I, C = Param_PulmonaryPeriphericalCirculation_C, V0 = Param_PulmonaryPeriphericalCirculation_V0, Vu0 = Param_PulmonaryPeriphericalCirculation_Vu0) annotation(
          Placement(visible = true, transformation(origin = {0.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Vessels.I1O1.Linear.CRL_LP PulmonaryVeins(R = Param_PulmonaryVeins_R, I = Param_PulmonaryVeins_I, C = Param_PulmonaryVeins_C, V0 = Param_PulmonaryVeins_V0, Vu0 = Param_PulmonaryVeins_Vu0) annotation(
          Placement(visible = true, transformation(origin = {50.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.AutoregulationCenters.Autoregulation_Center AutoregulationCenter(TauP = Param_AutoregulationCenter_TauP, TauZ = Param_AutoregulationCenter_TauZ, Pn = Param_AutoregulationCenter_Pn, ka = Param_AutoregulationCenter_ka, fmin = Param_AutoregulationCenter_fmin, fmax = Param_AutoregulationCenter_fmax, fesinf = Param_AutoregulationCenter_fesinf, fes0 = Param_AutoregulationCenter_fes0, fesmin = Param_AutoregulationCenter_fesmin, kes = Param_AutoregulationCenter_kes, fev0 = Param_AutoregulationCenter_fev0, fevinf = Param_AutoregulationCenter_fevinf, fcs0 = Param_AutoregulationCenter_fcs0, kev = Param_AutoregulationCenter_kev) annotation(
          Placement(visible = true, transformation(origin = {0.0, -30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -630)));
      equation
        connect(AutoregulationCenter.Active_fes, RightVentricle.Active_fes) annotation(
          Line(visible = true, points = {{-3.0, -22.0}, {-3.0, -13.0}, {-21.0, -13.0}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fes, SplanchnicPeripheralCirculation.Active_fes) annotation(
          Line(visible = true, origin = {-2.2398, -21.5986}, points = {{-0.76, -0.401}, {-0.76, 8.599}, {-11.994, 8.599}, {-11.994, -58.401}, {2.2398, -58.4014}, {2.24, -64.401}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fes, ExtraSplanchnicPeripheralCirculation.Active_fes) annotation(
          Line(visible = true, origin = {1.0, 0.0}, points = {{-4.0, -22.0}, {-4.0, -13.0}, {-16.0, -13.0}, {-16.0, -80.0}, {-1.0, -80.0}, {-1.0, -74.0}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fes, ExtraSplanchnicVeins.Active_fes) annotation(
          Line(visible = true, origin = {1.0, 0.0}, points = {{-4.0, -22.0}, {-4.0, -13.0}, {-16.0, -13.0}, {-16.0, -80.0}, {-71.0, -80.0}, {-71.0, -74.0}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fes, SplanchnicVeins.Active_fes) annotation(
          Line(visible = true, origin = {23.6786, 38.8776}, points = {{-26.679, -60.878}, {-26.679, -51.878}, {-38.679, -51.878}, {-38.679, -118.878}, {-93.679, -118.878}, {-93.679, -124.878}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fev, RightVentricle.Active_fev) annotation(
          Line(visible = true, points = {{3.0, -22.0}, {3.0, -7.0}, {-21.0, -7.0}}, color = {0, 255, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fev, LeftVentricle.Active_fev) annotation(
          Line(visible = true, points = {{3.0, -22.0}, {3.0, -7.0}, {21.0, -7.0}}, color = {0, 255, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fes, LeftVentricle.Active_fes) annotation(
          Line(visible = true, points = {{-3.0, -22.0}, {-3.0, -13.0}, {21.0, -13.0}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.PressureProbe, AorticValve.Outlet) annotation(
          Line(visible = true, points = {{0.0, -38.0}, {30.0, -38.0}}, color = {255, 0, 255}, thickness = 1));
        connect(ExtraSplanchnicVeins.Outlet, RightAtrium.Inlet) annotation(
          Line(visible = true, origin = {77.7551, 55.6165}, points = {{-155.755, -125.617}, {-167.755, -125.617}, {-167.755, -5.617}, {-107.755, -5.617}, {-107.755, -17.617}}, color = {0, 128, 255}, thickness = 2));
        connect(SplanchnicVeins.Outlet, RightAtrium.Inlet) annotation(
          Line(visible = true, origin = {-63.176, -2.1599}, points = {{-14.824, -87.84}, {-26.824, -87.84}, {-26.824, 52.16}, {33.176, 52.16}, {33.176, 40.16}}, color = {0, 128, 255}, thickness = 2));
        connect(ExtraSplanchnicPeripheralCirculation.Outlet, ExtraSplanchnicVeins.Inlet) annotation(
          Line(visible = true, origin = {-35.0, -70.0}, points = {{27.0, 0.0}, {-27.0, 0.0}}, color = {0, 128, 255}, thickness = 2));
        connect(SplanchnicPeripheralCirculation.Outlet, SplanchnicVeins.Inlet) annotation(
          Line(visible = true, origin = {-35.0, -90.0}, points = {{27.0, 0.0}, {-27.0, 0.0}}, color = {0, 128, 255}, thickness = 2));
        connect(SystemicArteries.Outlet, SplanchnicPeripheralCirculation.Inlet) annotation(
          Line(visible = true, points = {{62.0, -80.0}, {50.0, -80.0}, {50.0, -90.0}, {8.0, -90.0}}, color = {255, 0, 0}, thickness = 2));
        connect(SystemicArteries.Outlet, ExtraSplanchnicPeripheralCirculation.Inlet) annotation(
          Line(visible = true, origin = {50.0, -75.0}, points = {{12.0, -5.0}, {0.0, -5.0}, {0.0, 5.0}, {-42.0, 5.0}}, color = {255, 0, 0}, thickness = 2));
        connect(AorticValve.Outlet, SystemicArteries.Inlet) annotation(
          Line(visible = true, origin = {73.4354, -42.1173}, points = {{-43.435, 4.117}, {-43.435, -7.883}, {16.565, -7.883}, {16.565, -37.883}, {4.565, -37.883}}, color = {255, 0, 0}, thickness = 2));
        connect(PulmonaryVeins.Outlet, LeftAtrium.Inlet) annotation(
          Line(visible = true, points = {{58.0, 80.0}, {70.0, 80.0}, {70.0, 50.0}, {30.0, 50.0}, {30.0, 38.0}}, color = {255, 0, 0}, thickness = 2));
        connect(PulmonaryPeriphericalCirculation.Outlet, PulmonaryVeins.Inlet) annotation(
          Line(visible = true, origin = {25.0, 80.0}, points = {{-17.0, 0.0}, {17.0, 0.0}}, color = {255, 0, 0}, thickness = 2));
        connect(PulmonaryArteries.Outlet, PulmonaryPeriphericalCirculation.Inlet) annotation(
          Line(visible = true, origin = {-25.0, 80.0}, points = {{-17.0, 0.0}, {17.0, 0.0}}, color = {0, 128, 255}, thickness = 2));
        connect(PulmonaryValve.Outlet, PulmonaryArteries.Inlet) annotation(
          Line(visible = true, origin = {19.0, -76.0}, points = {{-49.0, 38.0}, {-49.0, 26.0}, {-89.0, 26.0}, {-89.0, 156.0}, {-77.0, 156.0}}, color = {0, 128, 255}, thickness = 2));
        connect(LeftVentricle.Outlet, AorticValve.Inlet) annotation(
          Line(visible = true, origin = {30.0, -20.6667}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
        connect(MitralicValve.Outlet, LeftVentricle.Inlet) annotation(
          Line(visible = true, origin = {30.0, -0.6667}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
        connect(LeftAtrium.Outlet, MitralicValve.Inlet) annotation(
          Line(visible = true, origin = {30.0, 19.3333}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
        connect(RightVentricle.Outlet, PulmonaryValve.Inlet) annotation(
          Line(visible = true, origin = {-30.0, -20.6667}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
        connect(TricuspidValve.Outlet, RightVentricle.Inlet) annotation(
          Line(visible = true, origin = {-30.0, -0.6667}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
        connect(RightAtrium.Outlet, TricuspidValve.Inlet) annotation(
          Line(visible = true, points = {{-30.0, 22.0}, {-30.0, 18.0}}));
        annotation(
          experiment(StopTime = 200, __Wolfram_SynchronizeWithRealTime = false),
          Diagram(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
      end Ursino1998Model;

      model ModelParametersNH
        parameter Real Param_LeftAtrium_C = 19.23;
        parameter Real Param_LeftAtrium_V0(unit = "ml") = 25;
        parameter Real Param_LeftAtrium_Vu0(unit = "ml") = 25;
        parameter Real Param_MitralicValve_R(unit = "mmHg.s2/ml") = 2.5 * 0.001;
        parameter Real Param_LeftVentricle_V0(unit = "ml") = 16.77;
        parameter Real Param_LeftVentricle_kR(unit = "s/ml") = 3.75 * 0.0001;
        parameter Real Param_LeftVentricle_Emax0 = 2.95;
        parameter Real Param_LeftVentricle_EmaxRef0 = 2.392;
        parameter Real Param_LeftVentricle_AGain_Emax = 0.475;
        parameter Real Param_LeftVentricle_ADelay_Emax = 2;
        parameter Real Param_LeftVentricle_ATau_Emax = 8;
        parameter Real Param_LeftVentricle_P0(unit = "mmHg") = 1.5;
        parameter Real Param_LeftVentricle_kE(unit = "1/ml") = 0.014;
        parameter Real Param_LeftVentricle_TSys0(unit = "s") = 0.5;
        parameter Real Param_LeftVentricle_kSys(unit = "s2") = 0.075;
        parameter Real Param_LeftVentricle_TRef0 = 0.58;
        parameter Real Param_LeftVentricle_AGain_Ts = -0.13;
        parameter Real Param_LeftVentricle_ADelay_Ts = 2;
        parameter Real Param_LeftVentricle_ATau_Ts = 2;
        parameter Real Param_LeftVentricle_AGain_Tv = 0.09;
        parameter Real Param_LeftVentricle_ADelay_Tv = 0.2;
        parameter Real Param_LeftVentricle_ATau_Tv = 1.5;
        parameter Real Param_AorticValve_kR = 3.75 * 0.0001;
        parameter Real Param_SystemicArteries_R(unit = "mmHg.s2/ml") = 0.06;
        parameter Real Param_SystemicArteries_I = 0.22 * 0.001;
        parameter Real Param_SystemicArteries_C(unit = "ml/mmHg") = 0.28;
        parameter Real Param_SystemicArteries_V0(unit = "ml") = 0;
        parameter Real Param_SystemicArteries_Vu0(unit = "ml") = 0;
        parameter Real Param_SplanchnicPeripheralCirculation_R0(unit = "mmHg.s2/ml") = 3.307;
        parameter Real Param_SplanchnicPeripheralCirculation_RRef0(unit = "mmHg.s2/ml") = 2.49;
        parameter Real Param_SplanchnicPeripheralCirculation_AGain = 0.695;
        parameter Real Param_SplanchnicPeripheralCirculation_ADelay = 2;
        parameter Real Param_SplanchnicPeripheralCirculation_ATau = 6;
        parameter Real Param_SplanchnicPeripheralCirculation_I(unit = "mmHg.s2/ml") = 0;
        parameter Real Param_SplanchnicPeripheralCirculation_C(unit = "ml/mmHg") = 2.05;
        parameter Real Param_SplanchnicPeripheralCirculation_V0(unit = "ml") = 274.4;
        parameter Real Param_SplanchnicPeripheralCirculation_Vu0(unit = "ml") = 274.4;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_R0(unit = "mmHg.s2/ml") = 1.407;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_RRef0(unit = "mmHg.s2/ml") = 0.78;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_AGain = 0.53;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_ADelay = 2;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_ATau = 6;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_I(unit = "mmHg.s2/ml") = 0;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_C(unit = "ml/mmHg") = 1.67;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_V0(unit = "ml") = 336.6;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_Vu0(unit = "ml") = 336.6;
        parameter Real Param_SplanchnicVeins_R(unit = "mmHg.s/ml") = 0.038;
        parameter Real Param_SplanchnicVeins_I(unit = "mmHg.s2/ml") = 0;
        parameter Real Param_SplanchnicVeins_C(unit = "ml/mmHg") = 61.11;
        parameter Real Param_SplanchnicVeins_V0 = 2963;
        parameter Real Param_SplanchnicVeins_Vu0 = 1121;
        parameter Real Param_SplanchnicVeins_VuRef0 = 1435.4;
        parameter Real Param_SplanchnicVeins_AGain = -265.4;
        parameter Real Param_SplanchnicVeins_ADelay = 5;
        parameter Real Param_SplanchnicVeins_ATau = 20;
        parameter Real Param_ExtraSplanchnicVeins_R(unit = "mmHg.s/ml") = 0.016;
        parameter Real Param_ExtraSplanchnicVeins_I(unit = "mmHg.s2/ml") = 0;
        parameter Real Param_ExtraSplanchnicVeins_C(unit = "ml/mmHg") = 50;
        parameter Real Param_ExtraSplanchnicVeins_V0(unit = "ml") = 1375;
        parameter Real Param_ExtraSplanchnicVeins_Vu0 = 1375;
        parameter Real Param_ExtraSplanchnicVeins_VuRef0(unit = "ml") = 1537;
        parameter Real Param_ExtraSplanchnicVeins_AGain = -132.5;
        parameter Real Param_ExtraSplanchnicVeins_ADelay = 5;
        parameter Real Param_ExtraSplanchnicVeins_ATau = 20;
        parameter Real Param_RightAtrium_C(unit = "ml/mmHg") = 31.25;
        parameter Real Param_RightAtrium_V0(unit = "ml") = 25;
        parameter Real Param_RightAtrium_Vu0(unit = "ml") = 25;
        parameter Real Param_TricuspidValve_R = 2.5 * 0.001;
        parameter Real Param_RightVentricle_V0(unit = "ml") = 40.8;
        parameter Real Param_RightVentricle_kR(unit = "s/ml") = 1.4 * 0.001;
        parameter Real Param_RightVentricle_Emax0 = 1.75;
        parameter Real Param_RightVentricle_EmaxRef0 = 1.412;
        parameter Real Param_RightVentricle_AGain_Emax = 0.282;
        parameter Real Param_RightVentricle_ADelay_Emax = 2;
        parameter Real Param_RightVentricle_ATau_Emax = 8;
        parameter Real Param_RightVentricle_P0(unit = "mmHg") = 1.5;
        parameter Real Param_RightVentricle_kE(unit = "1/ml") = 0.011;
        parameter Real Param_RightVentricle_TSys0(unit = "s") = 0.5;
        parameter Real Param_RightVentricle_kSys(unit = "s2") = 0.075;
        parameter Real Param_RightVentricle_TRef0 = 0.58;
        parameter Real Param_RightVentricle_AGain_Ts = -0.13;
        parameter Real Param_RightVentricle_ADelay_Ts = 2;
        parameter Real Param_RightVentricle_ATau_Ts = 2;
        parameter Real Param_RightVentricle_AGain_Tv = 0.09;
        parameter Real Param_RightVentricle_ADelay_Tv = 0.2;
        parameter Real Param_RightVentricle_ATau_Tv = 1.5;
        parameter Real Param_PulmonaryValve_kR = 1.4 * 0.001;
        parameter Real Param_PulmonaryArteries_R(unit = "mmHg.s/ml") = 0.023;
        parameter Real Param_PulmonaryArteries_I(unit = "mmHg.s2/ml") = 0.18 * 0.001;
        parameter Real Param_PulmonaryArteries_C(unit = "ml/mmHg") = 0.76;
        parameter Real Param_PulmonaryArteries_V0(unit = "ml") = 0;
        parameter Real Param_PulmonaryArteries_Vu0 = 0;
        parameter Real Param_PulmonaryPeriphericalCirculation_R(unit = "mmHg.s/ml") = 0.0894;
        parameter Real Param_PulmonaryPeriphericalCirculation_I(unit = "mmHg.s2/ml") = 0;
        parameter Real Param_PulmonaryPeriphericalCirculation_C(unit = "ml/mmHg") = 5.8;
        parameter Real Param_PulmonaryPeriphericalCirculation_V0(unit = "ml") = 123;
        parameter Real Param_PulmonaryPeriphericalCirculation_Vu0 = 123;
        parameter Real Param_PulmonaryVeins_R(unit = "mmHg.s/ml") = 0.0056;
        parameter Real Param_PulmonaryVeins_I(unit = "mmHg.s2/ml") = 0;
        parameter Real Param_PulmonaryVeins_C(unit = "ml/mmHg") = 25.37;
        parameter Real Param_PulmonaryVeins_V0(unit = "ml") = 120;
        parameter Real Param_PulmonaryVeins_Vu0(unit = "ml") = 120;
        parameter Real Param_AutoregulationCenter_TauP = 2.076;
        parameter Real Param_AutoregulationCenter_TauZ = 6.37;
        parameter Real Param_AutoregulationCenter_Pn = 92;
        parameter Real Param_AutoregulationCenter_ka = 11.758;
        parameter Real Param_AutoregulationCenter_fmin = 2.52;
        parameter Real Param_AutoregulationCenter_fmax = 47.78;
        parameter Real Param_AutoregulationCenter_fesinf = 2.1;
        parameter Real Param_AutoregulationCenter_fes0 = 16.11;
        parameter Real Param_AutoregulationCenter_fesmin = 2.66;
        parameter Real Param_AutoregulationCenter_kes = 0.0675;
        parameter Real Param_AutoregulationCenter_fev0 = 3.2;
        parameter Real Param_AutoregulationCenter_fevinf = 6.3;
        parameter Real Param_AutoregulationCenter_fcs0 = 25;
        parameter Real Param_AutoregulationCenter_kev = 7.06;
        parameter Real Param_LVAD_RPM(unit = "RPM") = 6000;
      end ModelParametersNH;

      model ModelParametersSHF_HTAP
        parameter Real Param_LeftAtrium_C = 19.23;
        parameter Real Param_LeftAtrium_V0(unit = "ml") = 25;
        parameter Real Param_LeftAtrium_Vu0(unit = "ml") = 25;
        parameter Real Param_MitralicValve_R(unit = "mmHg.s2/ml") = 2.5 * 0.001;
        parameter Real Param_LeftVentricle_V0(unit = "ml") = 16.77;
        parameter Real Param_LeftVentricle_kR(unit = "s/ml") = 3.75 * 0.0001;
        parameter Real Param_LeftVentricle_Emax0 = 0.2;
        parameter Real Param_LeftVentricle_EmaxRef0 = 0.2;
        parameter Real Param_LeftVentricle_AGain_Emax = 0.2;
        parameter Real Param_LeftVentricle_ADelay_Emax = 2;
        parameter Real Param_LeftVentricle_ATau_Emax = 8;
        parameter Real Param_LeftVentricle_P0(unit = "mmHg") = 1.5;
        parameter Real Param_LeftVentricle_kE(unit = "1/ml") = 0.011;
        parameter Real Param_LeftVentricle_TSys0(unit = "s") = 0.5;
        parameter Real Param_LeftVentricle_kSys(unit = "s2") = 0.075;
        parameter Real Param_LeftVentricle_TRef0 = 0.58;
        parameter Real Param_LeftVentricle_AGain_Ts = -0.13;
        parameter Real Param_LeftVentricle_ADelay_Ts = 2;
        parameter Real Param_LeftVentricle_ATau_Ts = 2;
        parameter Real Param_LeftVentricle_AGain_Tv = 0.09;
        parameter Real Param_LeftVentricle_ADelay_Tv = 0.2;
        parameter Real Param_LeftVentricle_ATau_Tv = 1.5;
        parameter Real Param_AorticValve_kR = 3.75 * 0.0001;
        parameter Real Param_SystemicArteries_R(unit = "mmHg.s2/ml") = 0.06;
        parameter Real Param_SystemicArteries_I = 0.22 * 0.001;
        parameter Real Param_SystemicArteries_C(unit = "ml/mmHg") = 0.28;
        parameter Real Param_SystemicArteries_V0(unit = "ml") = 0;
        parameter Real Param_SystemicArteries_Vu0(unit = "ml") = 0;
        parameter Real Param_SplanchnicPeripheralCirculation_R0(unit = "mmHg.s2/ml") = 3.307;
        parameter Real Param_SplanchnicPeripheralCirculation_RRef0(unit = "mmHg.s2/ml") = 2.49;
        parameter Real Param_SplanchnicPeripheralCirculation_AGain = 0.695;
        parameter Real Param_SplanchnicPeripheralCirculation_ADelay = 2;
        parameter Real Param_SplanchnicPeripheralCirculation_ATau = 6;
        parameter Real Param_SplanchnicPeripheralCirculation_I(unit = "mmHg.s2/ml") = 0;
        parameter Real Param_SplanchnicPeripheralCirculation_C(unit = "ml/mmHg") = 2.05;
        parameter Real Param_SplanchnicPeripheralCirculation_V0(unit = "ml") = 274.4;
        parameter Real Param_SplanchnicPeripheralCirculation_Vu0(unit = "ml") = 274.4;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_R0(unit = "mmHg.s2/ml") = 1.407;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_RRef0(unit = "mmHg.s2/ml") = 0.78;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_AGain = 0.53;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_ADelay = 2;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_ATau = 6;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_I(unit = "mmHg.s2/ml") = 0;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_C(unit = "ml/mmHg") = 1.67;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_V0(unit = "ml") = 336.6;
        parameter Real Param_ExtraSplanchnicPeripheralCirculation_Vu0(unit = "ml") = 336.6;
        parameter Real Param_SplanchnicVeins_R(unit = "mmHg.s/ml") = 0.038;
        parameter Real Param_SplanchnicVeins_I(unit = "mmHg.s2/ml") = 0;
        parameter Real Param_SplanchnicVeins_C(unit = "ml/mmHg") = 61.11;
        parameter Real Param_SplanchnicVeins_V0 = 2963;
        parameter Real Param_SplanchnicVeins_Vu0 = 1121;
        parameter Real Param_SplanchnicVeins_VuRef0 = 1435.4;
        parameter Real Param_SplanchnicVeins_AGain = -265.4;
        parameter Real Param_SplanchnicVeins_ADelay = 5;
        parameter Real Param_SplanchnicVeins_ATau = 20;
        parameter Real Param_ExtraSplanchnicVeins_R(unit = "mmHg.s/ml") = 0.016;
        parameter Real Param_ExtraSplanchnicVeins_I(unit = "mmHg.s2/ml") = 0;
        parameter Real Param_ExtraSplanchnicVeins_C(unit = "ml/mmHg") = 50;
        parameter Real Param_ExtraSplanchnicVeins_V0(unit = "ml") = 1375;
        parameter Real Param_ExtraSplanchnicVeins_Vu0 = 1375;
        parameter Real Param_ExtraSplanchnicVeins_VuRef0(unit = "ml") = 1537;
        parameter Real Param_ExtraSplanchnicVeins_AGain = -132.5;
        parameter Real Param_ExtraSplanchnicVeins_ADelay = 5;
        parameter Real Param_ExtraSplanchnicVeins_ATau = 20;
        parameter Real Param_RightAtrium_C(unit = "ml/mmHg") = 31.25;
        parameter Real Param_RightAtrium_V0(unit = "ml") = 25;
        parameter Real Param_RightAtrium_Vu0(unit = "ml") = 25;
        parameter Real Param_TricuspidValve_R = 2.5 * 0.001;
        parameter Real Param_RightVentricle_V0(unit = "ml") = 40.8;
        parameter Real Param_RightVentricle_kR(unit = "s/ml") = 1.4 * 0.001;
        parameter Real Param_RightVentricle_Emax0 = 1.75;
        parameter Real Param_RightVentricle_EmaxRef0 = 1.412;
        parameter Real Param_RightVentricle_AGain_Emax = 0.282;
        parameter Real Param_RightVentricle_ADelay_Emax = 2;
        parameter Real Param_RightVentricle_ATau_Emax = 8;
        parameter Real Param_RightVentricle_P0(unit = "mmHg") = 1.5;
        parameter Real Param_RightVentricle_kE(unit = "1/ml") = 0.011;
        parameter Real Param_RightVentricle_TSys0(unit = "s") = 0.5;
        parameter Real Param_RightVentricle_kSys(unit = "s2") = 0.075;
        parameter Real Param_RightVentricle_TRef0 = 0.58;
        parameter Real Param_RightVentricle_AGain_Ts = -0.13;
        parameter Real Param_RightVentricle_ADelay_Ts = 2;
        parameter Real Param_RightVentricle_ATau_Ts = 2;
        parameter Real Param_RightVentricle_AGain_Tv = 0.09;
        parameter Real Param_RightVentricle_ADelay_Tv = 0.2;
        parameter Real Param_RightVentricle_ATau_Tv = 1.5;
        parameter Real Param_PulmonaryValve_kR = 1.4 * 0.001;
        parameter Real Param_PulmonaryArteries_R(unit = "mmHg.s/ml") = 0.06;
        parameter Real Param_PulmonaryArteries_I(unit = "mmHg.s2/ml") = 0.18 * 0.001;
        parameter Real Param_PulmonaryArteries_C(unit = "ml/mmHg") = 0.76;
        parameter Real Param_PulmonaryArteries_V0(unit = "ml") = 0;
        parameter Real Param_PulmonaryArteries_Vu0 = 0;
        parameter Real Param_PulmonaryPeriphericalCirculation_R(unit = "mmHg.s/ml") = 0.4;
        parameter Real Param_PulmonaryPeriphericalCirculation_I(unit = "mmHg.s2/ml") = 0;
        parameter Real Param_PulmonaryPeriphericalCirculation_C(unit = "ml/mmHg") = 5.8;
        parameter Real Param_PulmonaryPeriphericalCirculation_V0(unit = "ml") = 123;
        parameter Real Param_PulmonaryPeriphericalCirculation_Vu0 = 123;
        parameter Real Param_PulmonaryVeins_R(unit = "mmHg.s/ml") = 0.0056;
        parameter Real Param_PulmonaryVeins_I(unit = "mmHg.s2/ml") = 0;
        parameter Real Param_PulmonaryVeins_C(unit = "ml/mmHg") = 25.37;
        parameter Real Param_PulmonaryVeins_V0(unit = "ml") = 120;
        parameter Real Param_PulmonaryVeins_Vu0(unit = "ml") = 120;
        parameter Real Param_AutoregulationCenter_TauP = 2.076;
        parameter Real Param_AutoregulationCenter_TauZ = 6.37;
        parameter Real Param_AutoregulationCenter_Pn = 92;
        parameter Real Param_AutoregulationCenter_ka = 11.758;
        parameter Real Param_AutoregulationCenter_fmin = 2.52;
        parameter Real Param_AutoregulationCenter_fmax = 47.78;
        parameter Real Param_AutoregulationCenter_fesinf = 2.1;
        parameter Real Param_AutoregulationCenter_fes0 = 16.11;
        parameter Real Param_AutoregulationCenter_fesmin = 2.66;
        parameter Real Param_AutoregulationCenter_kes = 0.0675;
        parameter Real Param_AutoregulationCenter_fev0 = 3.2;
        parameter Real Param_AutoregulationCenter_fevinf = 6.3;
        parameter Real Param_AutoregulationCenter_fcs0 = 25;
        parameter Real Param_AutoregulationCenter_kev = 7.06;
        parameter Real Param_LVAD_RPM(unit = "RPM") = 8000;
      end ModelParametersSHF_HTAP;

      model Ursino1998ModelECMO
        import Mathcard.Library.*;
        extends ModelParametersMHF10000;
        Devices.VAD LVAD(RPM = Param_LVAD_RPM) annotation(
          Placement(visible = true, transformation(origin = {80.0, -20.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Heart.AutoregulatingChambers.Atrium LeftAtrium(C = Param_LeftAtrium_C, V0 = Param_LeftAtrium_V0, Vu0 = Param_LeftAtrium_Vu0) annotation(
          Placement(visible = true, transformation(origin = {30.0, 30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Heart.Valves.Valve_AV MitralicValve(R = Param_MitralicValve_R) annotation(
          Placement(visible = true, transformation(origin = {30, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Heart.AutoregulatingChambers.Ventricle LeftVentricle(V0 = Param_LeftVentricle_V0, kR = Param_LeftVentricle_kR, Emax0 = Param_LeftVentricle_Emax0, EmaxRef0 = Param_LeftVentricle_EmaxRef0, AGain_Emax = Param_LeftVentricle_AGain_Emax, ADelay_Emax = Param_LeftVentricle_ADelay_Emax, ATau_Emax = Param_LeftVentricle_ATau_Emax, P0 = Param_LeftVentricle_P0, kE = Param_LeftVentricle_kE, TSys0 = Param_LeftVentricle_TSys0, kSys = Param_LeftVentricle_kSys, TRef0 = Param_LeftVentricle_TRef0, AGain_Ts = Param_LeftVentricle_AGain_Ts, ADelay_Ts = Param_LeftVentricle_ADelay_Ts, ATau_Ts = Param_LeftVentricle_ATau_Ts, AGain_Tv = Param_LeftVentricle_AGain_Tv, ADelay_Tv = Param_LeftVentricle_ADelay_Tv, ATau_Tv = Param_LeftVentricle_ATau_Tv) annotation(
          Placement(visible = true, transformation(origin = {30.0, -10.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
        Heart.Valves.Valve_VC AorticValve(kR = Param_AorticValve_kR) annotation(
          Placement(visible = true, transformation(origin = {30.0, -30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Vessels.I1O1.Linear.CRL_LP SystemicArteries(R = Param_SystemicArteries_R, I = Param_SystemicArteries_I, C = Param_SystemicArteries_C, V0 = Param_SystemicArteries_V0, Vu0 = Param_SystemicArteries_Vu0) annotation(
          Placement(visible = true, transformation(origin = {70.0, -80.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
        Vessels.I1O1.Autoregulating.CRL_LP_AR SplanchnicPeripheralCirculation(R0 = Param_SplanchnicPeripheralCirculation_R0, RRef0 = Param_SplanchnicPeripheralCirculation_RRef0, AGain = Param_SplanchnicPeripheralCirculation_AGain, ADelay = Param_SplanchnicPeripheralCirculation_ADelay, ATau = Param_SplanchnicPeripheralCirculation_ATau, I = Param_SplanchnicPeripheralCirculation_I, C = Param_SplanchnicPeripheralCirculation_C, V0 = Param_SplanchnicPeripheralCirculation_V0, Vu0 = Param_SplanchnicPeripheralCirculation_Vu0) annotation(
          Placement(visible = true, transformation(origin = {0.0, -90.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
        Vessels.I1O1.Autoregulating.CRL_LP_AR ExtraSplanchnicPeripheralCirculation(R0 = Param_ExtraSplanchnicPeripheralCirculation_R0, RRef0 = Param_ExtraSplanchnicPeripheralCirculation_RRef0, AGain = Param_ExtraSplanchnicPeripheralCirculation_AGain, ADelay = Param_ExtraSplanchnicPeripheralCirculation_ADelay, ATau = Param_ExtraSplanchnicPeripheralCirculation_ATau, I = Param_ExtraSplanchnicPeripheralCirculation_I, C = Param_ExtraSplanchnicPeripheralCirculation_C, V0 = Param_ExtraSplanchnicPeripheralCirculation_V0, Vu0 = Param_ExtraSplanchnicPeripheralCirculation_Vu0) annotation(
          Placement(visible = true, transformation(origin = {0.0, -70.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -180)));
        Vessels.I1O1.Autoregulating.CRL_LP_AVu SplanchnicVeins(R = Param_SplanchnicVeins_R, I = Param_SplanchnicVeins_I, C = Param_SplanchnicVeins_C, V0 = Param_SplanchnicVeins_V0, Vu0 = Param_SplanchnicVeins_Vu0, VuRef0 = Param_SplanchnicVeins_VuRef0, AGain = Param_SplanchnicVeins_AGain, ADelay = Param_SplanchnicVeins_ADelay, ATau = Param_SplanchnicVeins_ATau) annotation(
          Placement(visible = true, transformation(origin = {-70.0, -90.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
        Vessels.I1O1.Autoregulating.CRL_LP_AVu ExtraSplanchnicVeins(R = Param_ExtraSplanchnicVeins_R, I = Param_ExtraSplanchnicVeins_I, C = Param_ExtraSplanchnicVeins_C, V0 = Param_ExtraSplanchnicVeins_V0, Vu0 = Param_ExtraSplanchnicVeins_Vu0, VuRef0 = Param_ExtraSplanchnicVeins_VuRef0, AGain = Param_ExtraSplanchnicVeins_AGain, ADelay = Param_ExtraSplanchnicVeins_ADelay, ATau = Param_ExtraSplanchnicVeins_ATau) annotation(
          Placement(visible = true, transformation(origin = {-70.0, -70.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -180)));
        Heart.AutoregulatingChambers.Atrium RightAtrium(C = Param_RightAtrium_C, V0 = Param_RightAtrium_V0, Vu0 = Param_RightAtrium_Vu0) annotation(
          Placement(visible = true, transformation(origin = {-30.0, 30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Heart.Valves.Valve_AV TricuspidValve(R = Param_TricuspidValve_R) annotation(
          Placement(visible = true, transformation(origin = {-30.0, 10.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Heart.AutoregulatingChambers.Ventricle RightVentricle(V0 = Param_RightVentricle_V0, kR = Param_RightVentricle_kR, Emax0 = Param_RightVentricle_Emax0, EmaxRef0 = Param_RightVentricle_EmaxRef0, AGain_Emax = Param_RightVentricle_AGain_Emax, ADelay_Emax = Param_RightVentricle_ADelay_Emax, ATau_Emax = Param_RightVentricle_ATau_Emax, P0 = Param_RightVentricle_P0, kE = Param_RightVentricle_kE, TSys0 = Param_RightVentricle_TSys0, kSys = Param_RightVentricle_kSys, TRef0 = Param_RightVentricle_TRef0, AGain_Ts = Param_RightVentricle_AGain_Ts, ADelay_Ts = Param_RightVentricle_ADelay_Ts, ATau_Ts = Param_RightVentricle_ATau_Ts, AGain_Tv = Param_RightVentricle_AGain_Tv, ADelay_Tv = Param_RightVentricle_ADelay_Tv, ATau_Tv = Param_RightVentricle_ATau_Tv) annotation(
          Placement(visible = true, transformation(origin = {-30.0, -10.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Heart.Valves.Valve_VC PulmonaryValve(kR = Param_PulmonaryValve_kR) annotation(
          Placement(visible = true, transformation(origin = {-30.0, -30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Vessels.I1O1.Linear.CRL_LP PulmonaryArteries(R = Param_PulmonaryArteries_R, I = Param_PulmonaryArteries_I, C = Param_PulmonaryArteries_C, V0 = Param_PulmonaryArteries_V0, Vu0 = Param_PulmonaryArteries_Vu0) annotation(
          Placement(visible = true, transformation(origin = {-50.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -360)));
        Vessels.I1O1.Linear.CRL_LP PulmonaryPeriphericalCirculation(R = Param_PulmonaryPeriphericalCirculation_R, I = Param_PulmonaryPeriphericalCirculation_I, C = Param_PulmonaryPeriphericalCirculation_C, V0 = Param_PulmonaryPeriphericalCirculation_V0, Vu0 = Param_PulmonaryPeriphericalCirculation_Vu0) annotation(
          Placement(visible = true, transformation(origin = {0.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Vessels.I1O1.Linear.CRL_LP PulmonaryVeins(R = Param_PulmonaryVeins_R, I = Param_PulmonaryVeins_I, C = Param_PulmonaryVeins_C, V0 = Param_PulmonaryVeins_V0, Vu0 = Param_PulmonaryVeins_Vu0) annotation(
          Placement(visible = true, transformation(origin = {50.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        AutoregulationCenters.Autoregulation_Center AutoregulationCenter(TauP = Param_AutoregulationCenter_TauP, TauZ = Param_AutoregulationCenter_TauZ, Pn = Param_AutoregulationCenter_Pn, ka = Param_AutoregulationCenter_ka, fmin = Param_AutoregulationCenter_fmin, fmax = Param_AutoregulationCenter_fmax, fesinf = Param_AutoregulationCenter_fesinf, fes0 = Param_AutoregulationCenter_fes0, fesmin = Param_AutoregulationCenter_fesmin, kes = Param_AutoregulationCenter_kes, fev0 = Param_AutoregulationCenter_fev0, fevinf = Param_AutoregulationCenter_fevinf, fcs0 = Param_AutoregulationCenter_fcs0, kev = Param_AutoregulationCenter_kev) annotation(
          Placement(visible = true, transformation(origin = {0.0, -30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -630)));
        Devices.ECMO ECMO1(RPM = Param_LVAD_RPM) annotation(
          Placement(visible = true, transformation(origin = {92.514, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(LVAD.Inlet, MitralicValve.Outlet) annotation(
          Line(visible = true, points = {{80, -12}, {80, 2}, {30, 2}}, color = {255, 0, 0}, thickness = 2));
        connect(LVAD.Outlet, SystemicArteries.Inlet) annotation(
          Line(visible = true, points = {{80.0, -28.0}, {80.0, -50.0}, {90.0, -50.0}, {90.0, -80.0}, {78.0, -80.0}}, color = {255, 0, 0}, thickness = 2));
        connect(ECMO1.Inlet, TricuspidValve.Outlet) annotation(
          Line(visible = true, points = {{92.514, 28}, {80, 28}, {30, 2}}, color = {255, 0, 0}, thickness = 2));
        connect(ECMO1.Outlet, SystemicArteries.Inlet) annotation(
          Line(visible = true, points = {{92.514, 12}, {80, -50}, {90, -50}, {90, -80}, {78, -80}}, color = {255, 0, 0}, thickness = 2));
        connect(AutoregulationCenter.Active_fes, RightVentricle.Active_fes) annotation(
          Line(visible = true, points = {{-3.0, -22.0}, {-3.0, -13.0}, {-21.0, -13.0}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fes, SplanchnicPeripheralCirculation.Active_fes) annotation(
          Line(visible = true, origin = {-2.2398, -21.5986}, points = {{-0.76, -0.401}, {-0.76, 8.599}, {-11.994, 8.599}, {-11.994, -58.401}, {2.2398, -58.4014}, {2.24, -64.401}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fes, ExtraSplanchnicPeripheralCirculation.Active_fes) annotation(
          Line(visible = true, origin = {1.0, 0.0}, points = {{-4.0, -22.0}, {-4.0, -13.0}, {-16.0, -13.0}, {-16.0, -80.0}, {-1.0, -80.0}, {-1.0, -74.0}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fes, ExtraSplanchnicVeins.Active_fes) annotation(
          Line(visible = true, origin = {1.0, 0.0}, points = {{-4.0, -22.0}, {-4.0, -13.0}, {-16.0, -13.0}, {-16.0, -80.0}, {-71.0, -80.0}, {-71.0, -74.0}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fes, SplanchnicVeins.Active_fes) annotation(
          Line(visible = true, origin = {23.6786, 38.8776}, points = {{-26.679, -60.878}, {-26.679, -51.878}, {-38.679, -51.878}, {-38.679, -118.878}, {-93.679, -118.878}, {-93.679, -124.878}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fev, RightVentricle.Active_fev) annotation(
          Line(visible = true, points = {{3.0, -22.0}, {3.0, -7.0}, {-21.0, -7.0}}, color = {0, 255, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fev, LeftVentricle.Active_fev) annotation(
          Line(visible = true, points = {{3.0, -22.0}, {3.0, -7.0}, {21.0, -7.0}}, color = {0, 255, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fes, LeftVentricle.Active_fes) annotation(
          Line(visible = true, points = {{-3.0, -22.0}, {-3.0, -13.0}, {21.0, -13.0}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.PressureProbe, AorticValve.Outlet) annotation(
          Line(visible = true, points = {{0.0, -38.0}, {30.0, -38.0}}, color = {255, 0, 255}, thickness = 1));
        connect(ExtraSplanchnicVeins.Outlet, RightAtrium.Inlet) annotation(
          Line(visible = true, origin = {77.7551, 55.6165}, points = {{-155.755, -125.617}, {-167.755, -125.617}, {-167.755, -5.617}, {-107.755, -5.617}, {-107.755, -17.617}}, color = {0, 128, 255}, thickness = 2));
        connect(SplanchnicVeins.Outlet, RightAtrium.Inlet) annotation(
          Line(visible = true, origin = {-63.176, -2.1599}, points = {{-14.824, -87.84}, {-26.824, -87.84}, {-26.824, 52.16}, {33.176, 52.16}, {33.176, 40.16}}, color = {0, 128, 255}, thickness = 2));
        connect(ExtraSplanchnicPeripheralCirculation.Outlet, ExtraSplanchnicVeins.Inlet) annotation(
          Line(visible = true, origin = {-35.0, -70.0}, points = {{27.0, 0.0}, {-27.0, 0.0}}, color = {0, 128, 255}, thickness = 2));
        connect(SplanchnicPeripheralCirculation.Outlet, SplanchnicVeins.Inlet) annotation(
          Line(visible = true, origin = {-35.0, -90.0}, points = {{27.0, 0.0}, {-27.0, 0.0}}, color = {0, 128, 255}, thickness = 2));
        connect(SystemicArteries.Outlet, SplanchnicPeripheralCirculation.Inlet) annotation(
          Line(visible = true, points = {{62.0, -80.0}, {50.0, -80.0}, {50.0, -90.0}, {8.0, -90.0}}, color = {255, 0, 0}, thickness = 2));
        connect(SystemicArteries.Outlet, ExtraSplanchnicPeripheralCirculation.Inlet) annotation(
          Line(visible = true, origin = {50.0, -75.0}, points = {{12.0, -5.0}, {0.0, -5.0}, {0.0, 5.0}, {-42.0, 5.0}}, color = {255, 0, 0}, thickness = 2));
        connect(AorticValve.Outlet, SystemicArteries.Inlet) annotation(
          Line(visible = true, origin = {73.4354, -42.1173}, points = {{-43.435, 4.117}, {-43.435, -7.883}, {16.565, -7.883}, {16.565, -37.883}, {4.565, -37.883}}, color = {255, 0, 0}, thickness = 2));
        connect(PulmonaryVeins.Outlet, LeftAtrium.Inlet) annotation(
          Line(visible = true, points = {{58.0, 80.0}, {70.0, 80.0}, {70.0, 50.0}, {30.0, 50.0}, {30.0, 38.0}}, color = {255, 0, 0}, thickness = 2));
        connect(PulmonaryPeriphericalCirculation.Outlet, PulmonaryVeins.Inlet) annotation(
          Line(visible = true, origin = {25.0, 80.0}, points = {{-17.0, 0.0}, {17.0, 0.0}}, color = {255, 0, 0}, thickness = 2));
        connect(PulmonaryArteries.Outlet, PulmonaryPeriphericalCirculation.Inlet) annotation(
          Line(visible = true, origin = {-25.0, 80.0}, points = {{-17.0, 0.0}, {17.0, 0.0}}, color = {0, 128, 255}, thickness = 2));
        connect(PulmonaryValve.Outlet, PulmonaryArteries.Inlet) annotation(
          Line(visible = true, origin = {19.0, -76.0}, points = {{-49.0, 38.0}, {-49.0, 26.0}, {-89.0, 26.0}, {-89.0, 156.0}, {-77.0, 156.0}}, color = {0, 128, 255}, thickness = 2));
        connect(LeftVentricle.Outlet, AorticValve.Inlet) annotation(
          Line(visible = true, origin = {30.0, -20.6667}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
        connect(MitralicValve.Outlet, LeftVentricle.Inlet) annotation(
          Line(visible = true, origin = {30, -0.667}, points = {{-0, 2.667}, {-0, -1.333}, {0, -1.333}}));
        connect(LeftAtrium.Outlet, MitralicValve.Inlet) annotation(
          Line(visible = true, origin = {30, 19.333}, points = {{0, 2.667}, {-0, -1.333}, {-0, -1.333}}));
        connect(RightVentricle.Outlet, PulmonaryValve.Inlet) annotation(
          Line(visible = true, origin = {-30.0, -20.6667}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
        connect(TricuspidValve.Outlet, RightVentricle.Inlet) annotation(
          Line(visible = true, origin = {-30.0, -0.6667}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
        connect(RightAtrium.Outlet, TricuspidValve.Inlet) annotation(
          Line(visible = true, points = {{-30.0, 22.0}, {-30.0, 18.0}}));
        annotation(
          Diagram(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
      end Ursino1998ModelECMO;

      package HMIII
        model Ursino1998Model_VAD2 "HMIII VAD model"
          // EXTEND THE SUPERCLASS
          // CHOOSE THE HEART FAILURE LEVEL (MHF / SHF)
          //extends Mathcard.Applications.Ursino1998.HMIII.ModelParametersMHF;
          extends Mathcard.Applications.Ursino1998.ModelParametersNH;
          // IMPORT MATHCARD LIBRARY
          import Mathcard.Library.*;
          // ============================
          //
          // VAD DEVICE
          Mathcard.Library.Devices.VAD2 LVAD(RPM_c = Param_LVAD_RPM) annotation(
            Placement(visible = true, transformation(origin = {80.0, -20.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          // LEFT ATRIUM
          Mathcard.Library.Heart.AutoregulatingChambers.Atrium LeftAtrium(C = Param_LeftAtrium_C, V0 = Param_LeftAtrium_V0, Vu0 = Param_LeftAtrium_Vu0) annotation(
            Placement(visible = true, transformation(origin = {30.0, 30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          // MITRALIC VALVE
          Mathcard.Library.Heart.Valves.Valve_AV MitralicValve(R = Param_MitralicValve_R) annotation(
            Placement(visible = true, transformation(origin = {30.0, 10.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          // LEFT VENTRICLE
          Mathcard.Library.Heart.AutoregulatingChambers.Ventricle LeftVentricle(V0 = Param_LeftVentricle_V0, kR = Param_LeftVentricle_kR, Emax0 = Param_LeftVentricle_Emax0, EmaxRef0 = Param_LeftVentricle_EmaxRef0, AGain_Emax = Param_LeftVentricle_AGain_Emax, ADelay_Emax = Param_LeftVentricle_ADelay_Emax, ATau_Emax = Param_LeftVentricle_ATau_Emax, P0 = Param_LeftVentricle_P0, kE = Param_LeftVentricle_kE, TSys0 = Param_LeftVentricle_TSys0, kSys = Param_LeftVentricle_kSys, TRef0 = Param_LeftVentricle_TRef0, AGain_Ts = Param_LeftVentricle_AGain_Ts, ADelay_Ts = Param_LeftVentricle_ADelay_Ts, ATau_Ts = Param_LeftVentricle_ATau_Ts, AGain_Tv = Param_LeftVentricle_AGain_Tv, ADelay_Tv = Param_LeftVentricle_ADelay_Tv, ATau_Tv = Param_LeftVentricle_ATau_Tv) annotation(
            Placement(visible = true, transformation(origin = {30.0, -10.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
          // AORTIC VALVE
          Mathcard.Library.Heart.Valves.Valve_VC AorticValve(kR = Param_AorticValve_kR) annotation(
            Placement(visible = true, transformation(origin = {30.0, -30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          // SYSTEMIC ARTERIES
          Mathcard.Library.Vessels.I1O1.Linear.CRL_LP SystemicArteries(R = Param_SystemicArteries_R, I = Param_SystemicArteries_I, C = Param_SystemicArteries_C, V0 = Param_SystemicArteries_V0, Vu0 = Param_SystemicArteries_Vu0) annotation(
            Placement(visible = true, transformation(origin = {70.0, -80.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
          // SPLANCHNIC PERIPHERAL CIRCULATION
          Mathcard.Library.Vessels.I1O1.Autoregulating.CRL_LP_AR SplanchnicPeripheralCirculation(R0 = Param_SplanchnicPeripheralCirculation_R0, RRef0 = Param_SplanchnicPeripheralCirculation_RRef0, AGain = Param_SplanchnicPeripheralCirculation_AGain, ADelay = Param_SplanchnicPeripheralCirculation_ADelay, ATau = Param_SplanchnicPeripheralCirculation_ATau, I = Param_SplanchnicPeripheralCirculation_I, C = Param_SplanchnicPeripheralCirculation_C, V0 = Param_SplanchnicPeripheralCirculation_V0, Vu0 = Param_SplanchnicPeripheralCirculation_Vu0) annotation(
            Placement(visible = true, transformation(origin = {0.0, -90.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
          // EXTRASPLANCHNIC PERIPHERAL CIRCULATION
          Mathcard.Library.Vessels.I1O1.Autoregulating.CRL_LP_AR ExtraSplanchnicPeripheralCirculation(R0 = Param_ExtraSplanchnicPeripheralCirculation_R0, RRef0 = Param_ExtraSplanchnicPeripheralCirculation_RRef0, AGain = Param_ExtraSplanchnicPeripheralCirculation_AGain, ADelay = Param_ExtraSplanchnicPeripheralCirculation_ADelay, ATau = Param_ExtraSplanchnicPeripheralCirculation_ATau, I = Param_ExtraSplanchnicPeripheralCirculation_I, C = Param_ExtraSplanchnicPeripheralCirculation_C, V0 = Param_ExtraSplanchnicPeripheralCirculation_V0, Vu0 = Param_ExtraSplanchnicPeripheralCirculation_Vu0) annotation(
            Placement(visible = true, transformation(origin = {0.0, -70.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -180)));
          // SPLANCHNIC VEINS
          Mathcard.Library.Vessels.I1O1.Autoregulating.CRL_LP_AVu SplanchnicVeins(R = Param_SplanchnicVeins_R, I = Param_SplanchnicVeins_I, C = Param_SplanchnicVeins_C, V0 = Param_SplanchnicVeins_V0, Vu0 = Param_SplanchnicVeins_Vu0, VuRef0 = Param_SplanchnicVeins_VuRef0, AGain = Param_SplanchnicVeins_AGain, ADelay = Param_SplanchnicVeins_ADelay, ATau = Param_SplanchnicVeins_ATau) annotation(
            Placement(visible = true, transformation(origin = {-70.0, -90.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
          // EXTRASPLANCHNIC VEINS
          Mathcard.Library.Vessels.I1O1.Autoregulating.CRL_LP_AVu ExtraSplanchnicVeins(R = Param_ExtraSplanchnicVeins_R, I = Param_ExtraSplanchnicVeins_I, C = Param_ExtraSplanchnicVeins_C, V0 = Param_ExtraSplanchnicVeins_V0, Vu0 = Param_ExtraSplanchnicVeins_Vu0, VuRef0 = Param_ExtraSplanchnicVeins_VuRef0, AGain = Param_ExtraSplanchnicVeins_AGain, ADelay = Param_ExtraSplanchnicVeins_ADelay, ATau = Param_ExtraSplanchnicVeins_ATau) annotation(
            Placement(visible = true, transformation(origin = {-70.0, -70.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -180)));
          // RIGHT ATRIUM
          Mathcard.Library.Heart.AutoregulatingChambers.Atrium RightAtrium(C = Param_RightAtrium_C, V0 = Param_RightAtrium_V0, Vu0 = Param_RightAtrium_Vu0) annotation(
            Placement(visible = true, transformation(origin = {-30.0, 30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          // TRICUSPID VALVE
          Mathcard.Library.Heart.Valves.Valve_AV TricuspidValve(R = Param_TricuspidValve_R) annotation(
            Placement(visible = true, transformation(origin = {-30.0, 10.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          // RIGHT VENTRICLE
          Mathcard.Library.Heart.AutoregulatingChambers.Ventricle RightVentricle(V0 = Param_RightVentricle_V0, kR = Param_RightVentricle_kR, Emax0 = Param_RightVentricle_Emax0, EmaxRef0 = Param_RightVentricle_EmaxRef0, AGain_Emax = Param_RightVentricle_AGain_Emax, ADelay_Emax = Param_RightVentricle_ADelay_Emax, ATau_Emax = Param_RightVentricle_ATau_Emax, P0 = Param_RightVentricle_P0, kE = Param_RightVentricle_kE, TSys0 = Param_RightVentricle_TSys0, kSys = Param_RightVentricle_kSys, TRef0 = Param_RightVentricle_TRef0, AGain_Ts = Param_RightVentricle_AGain_Ts, ADelay_Ts = Param_RightVentricle_ADelay_Ts, ATau_Ts = Param_RightVentricle_ATau_Ts, AGain_Tv = Param_RightVentricle_AGain_Tv, ADelay_Tv = Param_RightVentricle_ADelay_Tv, ATau_Tv = Param_RightVentricle_ATau_Tv) annotation(
            Placement(visible = true, transformation(origin = {-30.0, -10.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          // PULMONARY VALVE
          Mathcard.Library.Heart.Valves.Valve_VC PulmonaryValve(kR = Param_PulmonaryValve_kR) annotation(
            Placement(visible = true, transformation(origin = {-30.0, -30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          // PULMONARY ARTERIES
          Mathcard.Library.Vessels.I1O1.Linear.CRL_LP PulmonaryArteries(R = Param_PulmonaryArteries_R, I = Param_PulmonaryArteries_I, C = Param_PulmonaryArteries_C, V0 = Param_PulmonaryArteries_V0, Vu0 = Param_PulmonaryArteries_Vu0) annotation(
            Placement(visible = true, transformation(origin = {-50.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -360)));
          // PULMONARY PERIPHERAL CIRCULATION
          Mathcard.Library.Vessels.I1O1.Linear.CRL_LP PulmonaryPeriphericalCirculation(R = Param_PulmonaryPeriphericalCirculation_R, I = Param_PulmonaryPeriphericalCirculation_I, C = Param_PulmonaryPeriphericalCirculation_C, V0 = Param_PulmonaryPeriphericalCirculation_V0, Vu0 = Param_PulmonaryPeriphericalCirculation_Vu0) annotation(
            Placement(visible = true, transformation(origin = {0.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          // PULMONARY VEINS
          Mathcard.Library.Vessels.I1O1.Linear.CRL_LP PulmonaryVeins(R = Param_PulmonaryVeins_R, I = Param_PulmonaryVeins_I, C = Param_PulmonaryVeins_C, V0 = Param_PulmonaryVeins_V0, Vu0 = Param_PulmonaryVeins_Vu0) annotation(
            Placement(visible = true, transformation(origin = {50.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
          // AUTOREGULATION CENTER
          Mathcard.Library.AutoregulationCenters.Autoregulation_Center AutoregulationCenter(TauP = Param_AutoregulationCenter_TauP, TauZ = Param_AutoregulationCenter_TauZ, Pn = Param_AutoregulationCenter_Pn, ka = Param_AutoregulationCenter_ka, fmin = Param_AutoregulationCenter_fmin, fmax = Param_AutoregulationCenter_fmax, fesinf = Param_AutoregulationCenter_fesinf, fes0 = Param_AutoregulationCenter_fes0, fesmin = Param_AutoregulationCenter_fesmin, kes = Param_AutoregulationCenter_kes, fev0 = Param_AutoregulationCenter_fev0, fevinf = Param_AutoregulationCenter_fevinf, fcs0 = Param_AutoregulationCenter_fcs0, kev = Param_AutoregulationCenter_kev) annotation(
            Placement(visible = true, transformation(origin = {0.0, -30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -630)));
        equation
          connect(LVAD.Inlet, MitralicValve.Outlet) annotation(
            Line(visible = true, points = {{80.0, -12.0}, {80.0, 2.0}, {30.0, 2.0}}, color = {255, 0, 0}, thickness = 2));
          connect(LVAD.Outlet, SystemicArteries.Inlet) annotation(
            Line(visible = true, points = {{80.0, -28.0}, {80.0, -50.0}, {90.0, -50.0}, {90.0, -80.0}, {78.0, -80.0}}, color = {255, 0, 0}, thickness = 2));
          connect(AutoregulationCenter.Active_fes, RightVentricle.Active_fes) annotation(
            Line(visible = true, points = {{-3.0, -22.0}, {-3.0, -13.0}, {-21.0, -13.0}}, color = {0, 128, 0}, thickness = 1));
          connect(AutoregulationCenter.Active_fes, SplanchnicPeripheralCirculation.Active_fes) annotation(
            Line(visible = true, origin = {-2.2398, -21.5986}, points = {{-0.76, -0.401}, {-0.76, 8.599}, {-11.994, 8.599}, {-11.994, -58.401}, {2.2398, -58.4014}, {2.24, -64.401}}, color = {0, 128, 0}, thickness = 1));
          connect(AutoregulationCenter.Active_fes, ExtraSplanchnicPeripheralCirculation.Active_fes) annotation(
            Line(visible = true, origin = {1.0, 0.0}, points = {{-4.0, -22.0}, {-4.0, -13.0}, {-16.0, -13.0}, {-16.0, -80.0}, {-1.0, -80.0}, {-1.0, -74.0}}, color = {0, 128, 0}, thickness = 1));
          connect(AutoregulationCenter.Active_fes, ExtraSplanchnicVeins.Active_fes) annotation(
            Line(visible = true, origin = {1.0, 0.0}, points = {{-4.0, -22.0}, {-4.0, -13.0}, {-16.0, -13.0}, {-16.0, -80.0}, {-71.0, -80.0}, {-71.0, -74.0}}, color = {0, 128, 0}, thickness = 1));
          connect(AutoregulationCenter.Active_fes, SplanchnicVeins.Active_fes) annotation(
            Line(visible = true, origin = {23.6786, 38.8776}, points = {{-26.679, -60.878}, {-26.679, -51.878}, {-38.679, -51.878}, {-38.679, -118.878}, {-93.679, -118.878}, {-93.679, -124.878}}, color = {0, 128, 0}, thickness = 1));
          connect(AutoregulationCenter.Active_fev, RightVentricle.Active_fev) annotation(
            Line(visible = true, points = {{3.0, -22.0}, {3.0, -7.0}, {-21.0, -7.0}}, color = {0, 255, 0}, thickness = 1));
          connect(AutoregulationCenter.Active_fev, LeftVentricle.Active_fev) annotation(
            Line(visible = true, points = {{3.0, -22.0}, {3.0, -7.0}, {21.0, -7.0}}, color = {0, 255, 0}, thickness = 1));
          connect(AutoregulationCenter.Active_fes, LeftVentricle.Active_fes) annotation(
            Line(visible = true, points = {{-3.0, -22.0}, {-3.0, -13.0}, {21.0, -13.0}}, color = {0, 128, 0}, thickness = 1));
          connect(AutoregulationCenter.PressureProbe, AorticValve.Outlet) annotation(
            Line(visible = true, points = {{0.0, -38.0}, {30.0, -38.0}}, color = {255, 0, 255}, thickness = 1));
          connect(ExtraSplanchnicVeins.Outlet, RightAtrium.Inlet) annotation(
            Line(visible = true, origin = {77.7551, 55.6165}, points = {{-155.755, -125.617}, {-167.755, -125.617}, {-167.755, -5.617}, {-107.755, -5.617}, {-107.755, -17.617}}, color = {0, 128, 255}, thickness = 2));
          connect(SplanchnicVeins.Outlet, RightAtrium.Inlet) annotation(
            Line(visible = true, origin = {-63.176, -2.1599}, points = {{-14.824, -87.84}, {-26.824, -87.84}, {-26.824, 52.16}, {33.176, 52.16}, {33.176, 40.16}}, color = {0, 128, 255}, thickness = 2));
          connect(ExtraSplanchnicPeripheralCirculation.Outlet, ExtraSplanchnicVeins.Inlet) annotation(
            Line(visible = true, origin = {-35.0, -70.0}, points = {{27.0, 0.0}, {-27.0, 0.0}}, color = {0, 128, 255}, thickness = 2));
          connect(SplanchnicPeripheralCirculation.Outlet, SplanchnicVeins.Inlet) annotation(
            Line(visible = true, origin = {-35.0, -90.0}, points = {{27.0, 0.0}, {-27.0, 0.0}}, color = {0, 128, 255}, thickness = 2));
          connect(SystemicArteries.Outlet, SplanchnicPeripheralCirculation.Inlet) annotation(
            Line(visible = true, points = {{62.0, -80.0}, {50.0, -80.0}, {50.0, -90.0}, {8.0, -90.0}}, color = {255, 0, 0}, thickness = 2));
          connect(SystemicArteries.Outlet, ExtraSplanchnicPeripheralCirculation.Inlet) annotation(
            Line(visible = true, origin = {50.0, -75.0}, points = {{12.0, -5.0}, {0.0, -5.0}, {0.0, 5.0}, {-42.0, 5.0}}, color = {255, 0, 0}, thickness = 2));
          connect(AorticValve.Outlet, SystemicArteries.Inlet) annotation(
            Line(visible = true, origin = {73.4354, -42.1173}, points = {{-43.435, 4.117}, {-43.435, -7.883}, {16.565, -7.883}, {16.565, -37.883}, {4.565, -37.883}}, color = {255, 0, 0}, thickness = 2));
          connect(PulmonaryVeins.Outlet, LeftAtrium.Inlet) annotation(
            Line(visible = true, points = {{58.0, 80.0}, {70.0, 80.0}, {70.0, 50.0}, {30.0, 50.0}, {30.0, 38.0}}, color = {255, 0, 0}, thickness = 2));
          connect(PulmonaryPeriphericalCirculation.Outlet, PulmonaryVeins.Inlet) annotation(
            Line(visible = true, origin = {25.0, 80.0}, points = {{-17.0, 0.0}, {17.0, 0.0}}, color = {255, 0, 0}, thickness = 2));
          connect(PulmonaryArteries.Outlet, PulmonaryPeriphericalCirculation.Inlet) annotation(
            Line(visible = true, origin = {-25.0, 80.0}, points = {{-17.0, 0.0}, {17.0, 0.0}}, color = {0, 128, 255}, thickness = 2));
          connect(PulmonaryValve.Outlet, PulmonaryArteries.Inlet) annotation(
            Line(visible = true, origin = {19.0, -76.0}, points = {{-49.0, 38.0}, {-49.0, 26.0}, {-89.0, 26.0}, {-89.0, 156.0}, {-77.0, 156.0}}, color = {0, 128, 255}, thickness = 2));
          connect(LeftVentricle.Outlet, AorticValve.Inlet) annotation(
            Line(visible = true, origin = {30.0, -20.6667}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
          connect(MitralicValve.Outlet, LeftVentricle.Inlet) annotation(
            Line(visible = true, origin = {30.0, -0.6667}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
          connect(LeftAtrium.Outlet, MitralicValve.Inlet) annotation(
            Line(visible = true, origin = {30.0, 19.3333}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
          connect(RightVentricle.Outlet, PulmonaryValve.Inlet) annotation(
            Line(visible = true, origin = {-30.0, -20.6667}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
          connect(TricuspidValve.Outlet, RightVentricle.Inlet) annotation(
            Line(visible = true, origin = {-30.0, -0.6667}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
          connect(RightAtrium.Outlet, TricuspidValve.Inlet) annotation(
            Line(visible = true, points = {{-30.0, 22.0}, {-30.0, 18.0}}));
          annotation(
            experiment(StopTime = 200),
            Diagram(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
        end Ursino1998Model_VAD2;

        partial model ModelParameters "Abstract class regrouping common parameters"
          /* PUMP RPMs */
          //parameter Real Param_LVAD_RPM = 3000;
          //parameter Real Param_LVAD_RPM = 4000;
          //parameter Real Param_LVAD_RPM = 5000;
          parameter Real Param_LVAD_RPM = 6000;
          //parameter Real Param_LVAD_RPM = 7000;
          //parameter Real Param_LVAD_RPM = 8000;
          //parameter Real Param_LVAD_RPM = 9000;
          //
          // ======================================
          //
          // Virtual parameters: related to heart failure level
          replaceable parameter Real Param_LeftVentricle_Emax0;
          replaceable parameter Real Param_LeftVentricle_EmaxRef0;
          replaceable parameter Real Param_LeftVentricle_kE(unit = "1/ml");
          //
          // ======================================
          //
          parameter Real Param_LeftAtrium_C = 19.23;
          parameter Real Param_LeftAtrium_V0(unit = "ml") = 25;
          parameter Real Param_LeftAtrium_Vu0(unit = "ml") = 25;
          parameter Real Param_MitralicValve_R(unit = "mmHg.s2/ml") = 2.5 * 0.001;
          parameter Real Param_LeftVentricle_V0(unit = "ml") = 16.77;
          parameter Real Param_LeftVentricle_kR(unit = "s/ml") = 3.75 * 0.0001;
          parameter Real Param_LeftVentricle_AGain_Emax = 0.2;
          parameter Real Param_LeftVentricle_ADelay_Emax = 2;
          parameter Real Param_LeftVentricle_ATau_Emax = 8;
          parameter Real Param_LeftVentricle_P0(unit = "mmHg") = 1.5;
          parameter Real Param_LeftVentricle_TSys0(unit = "s") = 0.5;
          parameter Real Param_LeftVentricle_kSys(unit = "s2") = 0.075;
          parameter Real Param_LeftVentricle_TRef0 = 0.58;
          parameter Real Param_LeftVentricle_AGain_Ts = -0.13;
          parameter Real Param_LeftVentricle_ADelay_Ts = 2;
          parameter Real Param_LeftVentricle_ATau_Ts = 2;
          parameter Real Param_LeftVentricle_AGain_Tv = 0.09;
          parameter Real Param_LeftVentricle_ADelay_Tv = 0.2;
          parameter Real Param_LeftVentricle_ATau_Tv = 1.5;
          parameter Real Param_AorticValve_kR = 3.75 * 0.0001;
          parameter Real Param_SystemicArteries_R(unit = "mmHg.s2/ml") = 0.06;
          parameter Real Param_SystemicArteries_I = 0.22 * 0.001;
          parameter Real Param_SystemicArteries_C(unit = "ml/mmHg") = 0.28;
          parameter Real Param_SystemicArteries_V0(unit = "ml") = 0;
          parameter Real Param_SystemicArteries_Vu0(unit = "ml") = 0;
          parameter Real Param_SplanchnicPeripheralCirculation_R0(unit = "mmHg.s2/ml") = 3.307;
          parameter Real Param_SplanchnicPeripheralCirculation_RRef0(unit = "mmHg.s2/ml") = 2.49;
          parameter Real Param_SplanchnicPeripheralCirculation_AGain = 0.695;
          parameter Real Param_SplanchnicPeripheralCirculation_ADelay = 2;
          parameter Real Param_SplanchnicPeripheralCirculation_ATau = 6;
          parameter Real Param_SplanchnicPeripheralCirculation_I(unit = "mmHg.s2/ml") = 0;
          parameter Real Param_SplanchnicPeripheralCirculation_C(unit = "ml/mmHg") = 2.05;
          parameter Real Param_SplanchnicPeripheralCirculation_V0(unit = "ml") = 274.4;
          parameter Real Param_SplanchnicPeripheralCirculation_Vu0(unit = "ml") = 274.4;
          parameter Real Param_ExtraSplanchnicPeripheralCirculation_R0(unit = "mmHg.s2/ml") = 1.407;
          parameter Real Param_ExtraSplanchnicPeripheralCirculation_RRef0(unit = "mmHg.s2/ml") = 0.78;
          parameter Real Param_ExtraSplanchnicPeripheralCirculation_AGain = 0.53;
          parameter Real Param_ExtraSplanchnicPeripheralCirculation_ADelay = 2;
          parameter Real Param_ExtraSplanchnicPeripheralCirculation_ATau = 6;
          parameter Real Param_ExtraSplanchnicPeripheralCirculation_I(unit = "mmHg.s2/ml") = 0;
          parameter Real Param_ExtraSplanchnicPeripheralCirculation_C(unit = "ml/mmHg") = 1.67;
          parameter Real Param_ExtraSplanchnicPeripheralCirculation_V0(unit = "ml") = 336.6;
          parameter Real Param_ExtraSplanchnicPeripheralCirculation_Vu0(unit = "ml") = 336.6;
          parameter Real Param_SplanchnicVeins_R(unit = "mmHg.s/ml") = 0.038;
          parameter Real Param_SplanchnicVeins_I(unit = "mmHg.s2/ml") = 0;
          parameter Real Param_SplanchnicVeins_C(unit = "ml/mmHg") = 61.11;
          parameter Real Param_SplanchnicVeins_V0 = 2963;
          parameter Real Param_SplanchnicVeins_Vu0 = 1121;
          parameter Real Param_SplanchnicVeins_VuRef0 = 1435.4;
          parameter Real Param_SplanchnicVeins_AGain = -265.4;
          parameter Real Param_SplanchnicVeins_ADelay = 5;
          parameter Real Param_SplanchnicVeins_ATau = 20;
          parameter Real Param_ExtraSplanchnicVeins_R(unit = "mmHg.s/ml") = 0.016;
          parameter Real Param_ExtraSplanchnicVeins_I(unit = "mmHg.s2/ml") = 0;
          parameter Real Param_ExtraSplanchnicVeins_C(unit = "ml/mmHg") = 50;
          parameter Real Param_ExtraSplanchnicVeins_V0(unit = "ml") = 1375;
          parameter Real Param_ExtraSplanchnicVeins_Vu0 = 1375;
          parameter Real Param_ExtraSplanchnicVeins_VuRef0(unit = "ml") = 1537;
          parameter Real Param_ExtraSplanchnicVeins_AGain = -132.5;
          parameter Real Param_ExtraSplanchnicVeins_ADelay = 5;
          parameter Real Param_ExtraSplanchnicVeins_ATau = 20;
          parameter Real Param_RightAtrium_C(unit = "ml/mmHg") = 31.25;
          parameter Real Param_RightAtrium_V0(unit = "ml") = 25;
          parameter Real Param_RightAtrium_Vu0(unit = "ml") = 25;
          parameter Real Param_TricuspidValve_R = 2.5 * 0.001;
          parameter Real Param_RightVentricle_V0(unit = "ml") = 40.8;
          parameter Real Param_RightVentricle_kR(unit = "s/ml") = 1.4 * 0.001;
          parameter Real Param_RightVentricle_Emax0 = 1.75;
          parameter Real Param_RightVentricle_EmaxRef0 = 1.412;
          parameter Real Param_RightVentricle_AGain_Emax = 0.282;
          parameter Real Param_RightVentricle_ADelay_Emax = 2;
          parameter Real Param_RightVentricle_ATau_Emax = 8;
          parameter Real Param_RightVentricle_P0(unit = "mmHg") = 1.5;
          parameter Real Param_RightVentricle_kE(unit = "1/ml") = 0.011;
          parameter Real Param_RightVentricle_TSys0(unit = "s") = 0.5;
          parameter Real Param_RightVentricle_kSys(unit = "s2") = 0.075;
          parameter Real Param_RightVentricle_TRef0 = 0.58;
          parameter Real Param_RightVentricle_AGain_Ts = -0.13;
          parameter Real Param_RightVentricle_ADelay_Ts = 2;
          parameter Real Param_RightVentricle_ATau_Ts = 2;
          parameter Real Param_RightVentricle_AGain_Tv = 0.09;
          parameter Real Param_RightVentricle_ADelay_Tv = 0.2;
          parameter Real Param_RightVentricle_ATau_Tv = 1.5;
          parameter Real Param_PulmonaryValve_kR = 1.4 * 0.001;
          parameter Real Param_PulmonaryArteries_R(unit = "mmHg.s/ml") = 0.023;
          parameter Real Param_PulmonaryArteries_I(unit = "mmHg.s2/ml") = 0.18 * 0.001;
          parameter Real Param_PulmonaryArteries_C(unit = "ml/mmHg") = 0.76;
          parameter Real Param_PulmonaryArteries_V0(unit = "ml") = 0;
          parameter Real Param_PulmonaryArteries_Vu0 = 0;
          parameter Real Param_PulmonaryPeriphericalCirculation_R(unit = "mmHg.s/ml") = 0.0894;
          parameter Real Param_PulmonaryPeriphericalCirculation_I(unit = "mmHg.s2/ml") = 0;
          parameter Real Param_PulmonaryPeriphericalCirculation_C(unit = "ml/mmHg") = 5.8;
          parameter Real Param_PulmonaryPeriphericalCirculation_V0(unit = "ml") = 123;
          parameter Real Param_PulmonaryPeriphericalCirculation_Vu0 = 123;
          parameter Real Param_PulmonaryVeins_R(unit = "mmHg.s/ml") = 0.0056;
          parameter Real Param_PulmonaryVeins_I(unit = "mmHg.s2/ml") = 0;
          parameter Real Param_PulmonaryVeins_C(unit = "ml/mmHg") = 25.37;
          parameter Real Param_PulmonaryVeins_V0(unit = "ml") = 120;
          parameter Real Param_PulmonaryVeins_Vu0(unit = "ml") = 120;
          parameter Real Param_AutoregulationCenter_TauP = 2.076;
          parameter Real Param_AutoregulationCenter_TauZ = 6.37;
          parameter Real Param_AutoregulationCenter_Pn = 92;
          parameter Real Param_AutoregulationCenter_ka = 11.758;
          parameter Real Param_AutoregulationCenter_fmin = 2.52;
          parameter Real Param_AutoregulationCenter_fmax = 47.78;
          parameter Real Param_AutoregulationCenter_fesinf = 2.1;
          parameter Real Param_AutoregulationCenter_fes0 = 16.11;
          parameter Real Param_AutoregulationCenter_fesmin = 2.66;
          parameter Real Param_AutoregulationCenter_kes = 0.0675;
          parameter Real Param_AutoregulationCenter_fev0 = 3.2;
          parameter Real Param_AutoregulationCenter_fevinf = 6.3;
          parameter Real Param_AutoregulationCenter_fcs0 = 25;
          parameter Real Param_AutoregulationCenter_kev = 7.06;
        end ModelParameters;

        model ModelParametersMHF "Medium-heart-failure-related parameters of left ventricle"
          // Inheritage
          extends ModelParameters;
          // Parameters
          redeclare parameter Real Param_LeftVentricle_Emax0 = 0.8;
          redeclare parameter Real Param_LeftVentricle_EmaxRef0 = 0.8;
          redeclare parameter Real Param_LeftVentricle_kE(unit = "1/ml") = 0.013;
        end ModelParametersMHF;

        model ModelParametersSHF "Severe-heart-failure-related parameters of left ventricle"
          // Inheritage
          extends ModelParameters;
          // Parameters
          redeclare parameter Real Param_LeftVentricle_Emax0 = 0.2;
          redeclare parameter Real Param_LeftVentricle_EmaxRef0 = 0.2;
          redeclare parameter Real Param_LeftVentricle_kE(unit = "1/ml") = 0.011;
        end ModelParametersSHF;
      end HMIII;

      model Ursino1998ModelVad
        //extends Mathcard.Applications.Ursino1998.ModelParametersNH;
        extends HMIII.ModelParametersMHF;
        import Mathcard.Library.*;
        Mathcard.Library.Devices.VAD LVAD(RPM = Param_LVAD_RPM) annotation(
          Placement(visible = true, transformation(origin = {80.0, -20.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Heart.AutoregulatingChambers.Atrium LeftAtrium(C = Param_LeftAtrium_C, V0 = Param_LeftAtrium_V0, Vu0 = Param_LeftAtrium_Vu0) annotation(
          Placement(visible = true, transformation(origin = {30.0, 30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Heart.Valves.Valve_AV MitralicValve(R = Param_MitralicValve_R) annotation(
          Placement(visible = true, transformation(origin = {30.0, 10.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Heart.AutoregulatingChambers.Ventricle LeftVentricle(V0 = Param_LeftVentricle_V0, kR = Param_LeftVentricle_kR, Emax0 = Param_LeftVentricle_Emax0, EmaxRef0 = Param_LeftVentricle_EmaxRef0, AGain_Emax = Param_LeftVentricle_AGain_Emax, ADelay_Emax = Param_LeftVentricle_ADelay_Emax, ATau_Emax = Param_LeftVentricle_ATau_Emax, P0 = Param_LeftVentricle_P0, kE = Param_LeftVentricle_kE, TSys0 = Param_LeftVentricle_TSys0, kSys = Param_LeftVentricle_kSys, TRef0 = Param_LeftVentricle_TRef0, AGain_Ts = Param_LeftVentricle_AGain_Ts, ADelay_Ts = Param_LeftVentricle_ADelay_Ts, ATau_Ts = Param_LeftVentricle_ATau_Ts, AGain_Tv = Param_LeftVentricle_AGain_Tv, ADelay_Tv = Param_LeftVentricle_ADelay_Tv, ATau_Tv = Param_LeftVentricle_ATau_Tv) annotation(
          Placement(visible = true, transformation(origin = {30.0, -10.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Heart.Valves.Valve_VC AorticValve(kR = Param_AorticValve_kR) annotation(
          Placement(visible = true, transformation(origin = {30.0, -30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Vessels.I1O1.Linear.CRL_LP SystemicArteries(R = Param_SystemicArteries_R, I = Param_SystemicArteries_I, C = Param_SystemicArteries_C, V0 = Param_SystemicArteries_V0, Vu0 = Param_SystemicArteries_Vu0) annotation(
          Placement(visible = true, transformation(origin = {70.0, -80.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Vessels.I1O1.Autoregulating.CRL_LP_AR SplanchnicPeripheralCirculation(R0 = Param_SplanchnicPeripheralCirculation_R0, RRef0 = Param_SplanchnicPeripheralCirculation_RRef0, AGain = Param_SplanchnicPeripheralCirculation_AGain, ADelay = Param_SplanchnicPeripheralCirculation_ADelay, ATau = Param_SplanchnicPeripheralCirculation_ATau, I = Param_SplanchnicPeripheralCirculation_I, C = Param_SplanchnicPeripheralCirculation_C, V0 = Param_SplanchnicPeripheralCirculation_V0, Vu0 = Param_SplanchnicPeripheralCirculation_Vu0) annotation(
          Placement(visible = true, transformation(origin = {0.0, -90.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Vessels.I1O1.Autoregulating.CRL_LP_AR ExtraSplanchnicPeripheralCirculation(R0 = Param_ExtraSplanchnicPeripheralCirculation_R0, RRef0 = Param_ExtraSplanchnicPeripheralCirculation_RRef0, AGain = Param_ExtraSplanchnicPeripheralCirculation_AGain, ADelay = Param_ExtraSplanchnicPeripheralCirculation_ADelay, ATau = Param_ExtraSplanchnicPeripheralCirculation_ATau, I = Param_ExtraSplanchnicPeripheralCirculation_I, C = Param_ExtraSplanchnicPeripheralCirculation_C, V0 = Param_ExtraSplanchnicPeripheralCirculation_V0, Vu0 = Param_ExtraSplanchnicPeripheralCirculation_Vu0) annotation(
          Placement(visible = true, transformation(origin = {0.0, -70.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -180)));
        Mathcard.Library.Vessels.I1O1.Autoregulating.CRL_LP_AVu SplanchnicVeins(R = Param_SplanchnicVeins_R, I = Param_SplanchnicVeins_I, C = Param_SplanchnicVeins_C, V0 = Param_SplanchnicVeins_V0, Vu0 = Param_SplanchnicVeins_Vu0, VuRef0 = Param_SplanchnicVeins_VuRef0, AGain = Param_SplanchnicVeins_AGain, ADelay = Param_SplanchnicVeins_ADelay, ATau = Param_SplanchnicVeins_ATau) annotation(
          Placement(visible = true, transformation(origin = {-70.0, -90.0}, extent = {{10.0, -10.0}, {-10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Vessels.I1O1.Autoregulating.CRL_LP_AVu ExtraSplanchnicVeins(R = Param_ExtraSplanchnicVeins_R, I = Param_ExtraSplanchnicVeins_I, C = Param_ExtraSplanchnicVeins_C, V0 = Param_ExtraSplanchnicVeins_V0, Vu0 = Param_ExtraSplanchnicVeins_Vu0, VuRef0 = Param_ExtraSplanchnicVeins_VuRef0, AGain = Param_ExtraSplanchnicVeins_AGain, ADelay = Param_ExtraSplanchnicVeins_ADelay, ATau = Param_ExtraSplanchnicVeins_ATau) annotation(
          Placement(visible = true, transformation(origin = {-70.0, -70.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -180)));
        Mathcard.Library.Heart.AutoregulatingChambers.Atrium RightAtrium(C = Param_RightAtrium_C, V0 = Param_RightAtrium_V0, Vu0 = Param_RightAtrium_Vu0) annotation(
          Placement(visible = true, transformation(origin = {-30.0, 30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Heart.Valves.Valve_AV TricuspidValve(R = Param_TricuspidValve_R) annotation(
          Placement(visible = true, transformation(origin = {-30.0, 10.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Heart.AutoregulatingChambers.Ventricle RightVentricle(V0 = Param_RightVentricle_V0, kR = Param_RightVentricle_kR, Emax0 = Param_RightVentricle_Emax0, EmaxRef0 = Param_RightVentricle_EmaxRef0, AGain_Emax = Param_RightVentricle_AGain_Emax, ADelay_Emax = Param_RightVentricle_ADelay_Emax, ATau_Emax = Param_RightVentricle_ATau_Emax, P0 = Param_RightVentricle_P0, kE = Param_RightVentricle_kE, TSys0 = Param_RightVentricle_TSys0, kSys = Param_RightVentricle_kSys, TRef0 = Param_RightVentricle_TRef0, AGain_Ts = Param_RightVentricle_AGain_Ts, ADelay_Ts = Param_RightVentricle_ADelay_Ts, ATau_Ts = Param_RightVentricle_ATau_Ts, AGain_Tv = Param_RightVentricle_AGain_Tv, ADelay_Tv = Param_RightVentricle_ADelay_Tv, ATau_Tv = Param_RightVentricle_ATau_Tv) annotation(
          Placement(visible = true, transformation(origin = {-30.0, -10.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Heart.Valves.Valve_VC PulmonaryValve(kR = Param_PulmonaryValve_kR) annotation(
          Placement(visible = true, transformation(origin = {-30.0, -30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Vessels.I1O1.Linear.CRL_LP PulmonaryArteries(R = Param_PulmonaryArteries_R, I = Param_PulmonaryArteries_I, C = Param_PulmonaryArteries_C, V0 = Param_PulmonaryArteries_V0, Vu0 = Param_PulmonaryArteries_Vu0) annotation(
          Placement(visible = true, transformation(origin = {-50.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -360)));
        Mathcard.Library.Vessels.I1O1.Linear.CRL_LP PulmonaryPeriphericalCirculation(R = Param_PulmonaryPeriphericalCirculation_R, I = Param_PulmonaryPeriphericalCirculation_I, C = Param_PulmonaryPeriphericalCirculation_C, V0 = Param_PulmonaryPeriphericalCirculation_V0, Vu0 = Param_PulmonaryPeriphericalCirculation_Vu0) annotation(
          Placement(visible = true, transformation(origin = {0.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.Vessels.I1O1.Linear.CRL_LP PulmonaryVeins(R = Param_PulmonaryVeins_R, I = Param_PulmonaryVeins_I, C = Param_PulmonaryVeins_C, V0 = Param_PulmonaryVeins_V0, Vu0 = Param_PulmonaryVeins_Vu0) annotation(
          Placement(visible = true, transformation(origin = {50.0, 80.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = 0)));
        Mathcard.Library.AutoregulationCenters.Autoregulation_Center AutoregulationCenter(TauP = Param_AutoregulationCenter_TauP, TauZ = Param_AutoregulationCenter_TauZ, Pn = Param_AutoregulationCenter_Pn, ka = Param_AutoregulationCenter_ka, fmin = Param_AutoregulationCenter_fmin, fmax = Param_AutoregulationCenter_fmax, fesinf = Param_AutoregulationCenter_fesinf, fes0 = Param_AutoregulationCenter_fes0, fesmin = Param_AutoregulationCenter_fesmin, kes = Param_AutoregulationCenter_kes, fev0 = Param_AutoregulationCenter_fev0, fevinf = Param_AutoregulationCenter_fevinf, fcs0 = Param_AutoregulationCenter_fcs0, kev = Param_AutoregulationCenter_kev) annotation(
          Placement(visible = true, transformation(origin = {0.0, -30.0}, extent = {{-10.0, -10.0}, {10.0, 10.0}}, rotation = -630)));
      equation
        connect(LVAD.Inlet, MitralicValve.Outlet) annotation(
          Line(visible = true, points = {{80.0, -12.0}, {80.0, 2.0}, {30.0, 2.0}}, color = {255, 0, 0}, thickness = 2));
        connect(LVAD.Outlet, SystemicArteries.Inlet) annotation(
          Line(visible = true, points = {{80.0, -28.0}, {80.0, -50.0}, {90.0, -50.0}, {90.0, -80.0}, {78.0, -80.0}}, color = {255, 0, 0}, thickness = 2));
        connect(AutoregulationCenter.Active_fes, RightVentricle.Active_fes) annotation(
          Line(visible = true, points = {{-3.0, -22.0}, {-3.0, -13.0}, {-21.0, -13.0}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fes, SplanchnicPeripheralCirculation.Active_fes) annotation(
          Line(visible = true, origin = {-2.2398, -21.5986}, points = {{-0.76, -0.401}, {-0.76, 8.599}, {-11.994, 8.599}, {-11.994, -58.401}, {2.2398, -58.4014}, {2.24, -64.401}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fes, ExtraSplanchnicPeripheralCirculation.Active_fes) annotation(
          Line(visible = true, origin = {1.0, 0.0}, points = {{-4.0, -22.0}, {-4.0, -13.0}, {-16.0, -13.0}, {-16.0, -80.0}, {-1.0, -80.0}, {-1.0, -74.0}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fes, ExtraSplanchnicVeins.Active_fes) annotation(
          Line(visible = true, origin = {1.0, 0.0}, points = {{-4.0, -22.0}, {-4.0, -13.0}, {-16.0, -13.0}, {-16.0, -80.0}, {-71.0, -80.0}, {-71.0, -74.0}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fes, SplanchnicVeins.Active_fes) annotation(
          Line(visible = true, origin = {23.6786, 38.8776}, points = {{-26.679, -60.878}, {-26.679, -51.878}, {-38.679, -51.878}, {-38.679, -118.878}, {-93.679, -118.878}, {-93.679, -124.878}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fev, RightVentricle.Active_fev) annotation(
          Line(visible = true, points = {{3.0, -22.0}, {3.0, -7.0}, {-21.0, -7.0}}, color = {0, 255, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fev, LeftVentricle.Active_fev) annotation(
          Line(visible = true, points = {{3.0, -22.0}, {3.0, -7.0}, {21.0, -7.0}}, color = {0, 255, 0}, thickness = 1));
        connect(AutoregulationCenter.Active_fes, LeftVentricle.Active_fes) annotation(
          Line(visible = true, points = {{-3.0, -22.0}, {-3.0, -13.0}, {21.0, -13.0}}, color = {0, 128, 0}, thickness = 1));
        connect(AutoregulationCenter.PressureProbe, AorticValve.Outlet) annotation(
          Line(visible = true, points = {{0.0, -38.0}, {30.0, -38.0}}, color = {255, 0, 255}, thickness = 1));
        connect(ExtraSplanchnicVeins.Outlet, RightAtrium.Inlet) annotation(
          Line(visible = true, origin = {77.7551, 55.6165}, points = {{-155.755, -125.617}, {-167.755, -125.617}, {-167.755, -5.617}, {-107.755, -5.617}, {-107.755, -17.617}}, color = {0, 128, 255}, thickness = 2));
        connect(SplanchnicVeins.Outlet, RightAtrium.Inlet) annotation(
          Line(visible = true, origin = {-63.176, -2.1599}, points = {{-14.824, -87.84}, {-26.824, -87.84}, {-26.824, 52.16}, {33.176, 52.16}, {33.176, 40.16}}, color = {0, 128, 255}, thickness = 2));
        connect(ExtraSplanchnicPeripheralCirculation.Outlet, ExtraSplanchnicVeins.Inlet) annotation(
          Line(visible = true, origin = {-35.0, -70.0}, points = {{27.0, 0.0}, {-27.0, 0.0}}, color = {0, 128, 255}, thickness = 2));
        connect(SplanchnicPeripheralCirculation.Outlet, SplanchnicVeins.Inlet) annotation(
          Line(visible = true, origin = {-35.0, -90.0}, points = {{27.0, 0.0}, {-27.0, 0.0}}, color = {0, 128, 255}, thickness = 2));
        connect(SystemicArteries.Outlet, SplanchnicPeripheralCirculation.Inlet) annotation(
          Line(visible = true, points = {{62.0, -80.0}, {50.0, -80.0}, {50.0, -90.0}, {8.0, -90.0}}, color = {255, 0, 0}, thickness = 2));
        connect(SystemicArteries.Outlet, ExtraSplanchnicPeripheralCirculation.Inlet) annotation(
          Line(visible = true, origin = {50.0, -75.0}, points = {{12.0, -5.0}, {0.0, -5.0}, {0.0, 5.0}, {-42.0, 5.0}}, color = {255, 0, 0}, thickness = 2));
        connect(AorticValve.Outlet, SystemicArteries.Inlet) annotation(
          Line(visible = true, origin = {73.4354, -42.1173}, points = {{-43.435, 4.117}, {-43.435, -7.883}, {16.565, -7.883}, {16.565, -37.883}, {4.565, -37.883}}, color = {255, 0, 0}, thickness = 2));
        connect(PulmonaryVeins.Outlet, LeftAtrium.Inlet) annotation(
          Line(visible = true, points = {{58.0, 80.0}, {70.0, 80.0}, {70.0, 50.0}, {30.0, 50.0}, {30.0, 38.0}}, color = {255, 0, 0}, thickness = 2));
        connect(PulmonaryPeriphericalCirculation.Outlet, PulmonaryVeins.Inlet) annotation(
          Line(visible = true, origin = {25.0, 80.0}, points = {{-17.0, 0.0}, {17.0, 0.0}}, color = {255, 0, 0}, thickness = 2));
        connect(PulmonaryArteries.Outlet, PulmonaryPeriphericalCirculation.Inlet) annotation(
          Line(visible = true, origin = {-25.0, 80.0}, points = {{-17.0, 0.0}, {17.0, 0.0}}, color = {0, 128, 255}, thickness = 2));
        connect(PulmonaryValve.Outlet, PulmonaryArteries.Inlet) annotation(
          Line(visible = true, origin = {19.0, -76.0}, points = {{-49.0, 38.0}, {-49.0, 26.0}, {-89.0, 26.0}, {-89.0, 156.0}, {-77.0, 156.0}}, color = {0, 128, 255}, thickness = 2));
        connect(LeftVentricle.Outlet, AorticValve.Inlet) annotation(
          Line(visible = true, origin = {30.0, -20.6667}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
        connect(MitralicValve.Outlet, LeftVentricle.Inlet) annotation(
          Line(visible = true, origin = {30.0, -0.6667}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
        connect(LeftAtrium.Outlet, MitralicValve.Inlet) annotation(
          Line(visible = true, origin = {30.0, 19.3333}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
        connect(RightVentricle.Outlet, PulmonaryValve.Inlet) annotation(
          Line(visible = true, origin = {-30.0, -20.6667}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
        connect(TricuspidValve.Outlet, RightVentricle.Inlet) annotation(
          Line(visible = true, origin = {-30.0, -0.6667}, points = {{0.0, 2.6667}, {0.0, -1.3333}, {0.0, -1.3333}}));
        connect(RightAtrium.Outlet, TricuspidValve.Inlet) annotation(
          Line(visible = true, points = {{-30.0, 22.0}, {-30.0, 18.0}}));
        annotation(
          experiment(StopTime = 20, Interval = 0.002),
          Diagram(coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
      end Ursino1998ModelVad;
      annotation(
        Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Ursino 1998", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
        Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
        Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
    end Ursino1998;
    annotation(
      Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Applications", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
      Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
  end Applications;
  annotation(
    Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {80, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{-100, 50}, {-80, 70}, {100, 70}, {80, 50}, {-100, 50}}), Polygon(visible = true, fillColor = {210, 210, 210}, fillPattern = FillPattern.Solid, points = {{100, 70}, {100, -80}, {80, -100}, {80, 50}, {100, 70}}), Text(visible = true, fillColor = {0, 0, 255}, extent = {{-85, -85}, {65, 35}}, textString = "Mathcard", fontName = "Arial"), Text(visible = true, fillColor = {255, 0, 0}, extent = {{-120, 73}, {120, 122}}, textString = "%name", fontName = "Arial")}),
    Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})),
    Diagram(coordinateSystem(extent = {{-148.5, -105.0}, {148.5, 105.0}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end Mathcard;