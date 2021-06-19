var Body = function (vertices, center, theta, mInv, moiInv, vLin, vAng) {
    this.vertices = vertices; // array of vertices relative to center
    this.center = center; // center of mass
    this.theta = theta;   // rotation angle
    this.mInv = mInv;     // inverse of the mass
    this.moiInv = moiInv; // inverse of the moment of inertia
    this.vLin = vLin;     // linear (translational) velocity
    this.vAng = vAng;     // angular (rotational) velocity
    this.forces = [];     // array of forces...
    this.forcePoints = [];// ... and the vertex index of force application.
    // undefined is center of mass.

    var rotMatTheta;      // used to avoid unnecessary rotation matrix computations
    var rotationMatrix;

    this.addForce = function (force, forcePoint) {
        this.forces.push(force);
        this.forcePoints.push(forcePoint);
    };
    this.getRotationMatrix = function () {
        // only recompute if theta has changed since the last call.
        if (this.theta !== rotMatTheta) {
            rotationMatrix = [[Math.cos(this.theta), -Math.sin(this.theta)],
            [Math.sin(this.theta), Math.cos(this.theta)]];
            rotMatTheta = this.theta;
        }
        return rotationMatrix;
    };
    this.getVerticesInWorldCoords = function () {
        var vertsAbsolute = [];
        var rotationMatrix = this.getRotationMatrix();
        for (var i = 0; i < this.vertices.length; i++) {
            vertsAbsolute.push(MV.VpV(this.center,
                MV.MxV(rotationMatrix, this.vertices[i])));
        }
        return vertsAbsolute;
    };

    this.getRotationMatrix();
};

var Joint = function (bodyA, vertexA, bodyB, vertexB) {
    this.bodyA = bodyA;
    this.vertexA = vertexA;
    this.bodyB = bodyB;
    this.vertexB = vertexB;
};


// world
// helpers
// step
// vel updater
// constraint solver

var World = function () {
    this.dt = 1 / 240;
    this.nIterations = 4;
    this.beta = 0.2;
    this.bodies = [];
    this.t = 0;
    this.joints = [];
    this.step = function () {
        /***** 1. integrate forces/torques and compute tentative velocities *****/
        for (var i = 0; i < this.bodies.length; i++) {
            var mInv = this.bodies[i].mInv;
            if (mInv == 0) continue;
            var moiInv = this.bodies[i].moiInv;
            var rotationMatrix = this.bodies[i].getRotationMatrix();
            // each force independently leads to a change in linear and angular
            // velocity; i.e., add up all these changes by iterating over forces
            for (var j = 0; j < this.bodies[i].forces.length; j++) {
                var force = this.bodies[i].forces[j];
                var forcePoint = this.bodies[i].forcePoints[j];
                // linear motion is simply the integrated force, divided by the mass
                this.bodies[i].vLin = MV.VpV(this.bodies[i].vLin, MV.SxV(this.dt * mInv, force));

                ////////////////////////////////////////////////////////////////////////////////////
                // should be split into linear & angular components for off center forces. below calculation isn't right i think. applies all newtowns above and doubles up down below. excess force is applied.
                ////////////////////////////////////////////////////////////////////////////////////

                // angular motion depends on the force application point as well via the torque
                if (forcePoint !== undefined) { // 'undefined' means center of mass
                    var torque = MV.cross2(MV.MxV(rotationMatrix, this.bodies[i].vertices[forcePoint]), force);
                    this.bodies[i].vAng += this.dt * moiInv * torque;
                }
            }
        }

        /***** 2. correct velocity errors *****/
        applyImpulses.apply(this);

        /***** 3. update positions *****/
        for (var i = 0; i < this.bodies.length; i++) {
            this.bodies[i].center = MV.VpV(this.bodies[i].center, MV.SxV(this.dt, this.bodies[i].vLin));
            this.bodies[i].theta += this.dt * this.bodies[i].vAng;
        }
    };

    this.addBody = function (body) {
        this.bodies.push(body);
    };

    this.addJoint = function (joint) {
        this.joints.push(joint);
    };

    function applyImpulses() {
        var bias = [];
        // precompute MInv and bias for each contact - they don't change
        // across iterations
        var MInv = [];
        var bias = [];
        var J = [];

        for (var i = 0; i < this.joints.length; i++) {
            // compute J, C, and bias
            var joint = this.joints[i];
            MInv[i] = [
                joint.bodyA.mInv, joint.bodyA.mInv, joint.bodyA.moiInv,
                joint.bodyB.mInv, joint.bodyB.mInv, joint.bodyB.moiInv,
            ]

            var pA = joint.bodyA.getVerticesInWorldCoords()[joint.vertexA] // contact point on body A (world coords)
            var pB = joint.bodyB.getVerticesInWorldCoords()[joint.vertexB] // contact point on body B (world coords)
            var cA = joint.bodyA.center // center of body A
            var cB = joint.bodyB.center // center of body B

            J[i] = MV.SxV(2,
                MV.VmV(pA, pB)
                .concat(MV.cross2(MV.VmV(pB, pA), MV.VmV(pA, cA)))
                .concat(MV.VmV(pB, pA))
                .concat(MV.cross2(MV.VmV(pA, pB), MV.VmV(pB, cB)))
            )

            C = MV.dot(MV.VmV(pA, pB), MV.VmV(pA, pB))
            bias[i] = this.beta / this.dt * C
        };

        for (var i = 0; i < this.nIterations; i++) {
            for (var j = 0; j < this.joints.length; j++) {
                // apply forces to both bodies in joint using MInv * lambda * J
                var bodyA = this.joints[j].bodyA
                var bodyB = this.joints[j].bodyB

                var v = bodyA.vLin
                        .concat(bodyA.vAng)
                        .concat(bodyB.vLin)
                        .concat(bodyB.vAng)

                var lambdaDenominator = MV.dot(J[j], MV.VxV(MInv[j], J[j]))
                if (Math.abs(lambdaDenominator) <= 1e-15) continue;
                var lambda = -(MV.dot(J[j], v) + bias[j]) / lambdaDenominator

                v = MV.VpV(
                    v,
                    MV.SxV(
                        lambda, MV.VxV(J[j], MInv[j])
                    )
                )

                bodyA.vLin = v.slice(0,2)
                bodyA.vAng = v[2]
                bodyB.vLin = v.slice(3,5)
                bodyB.vAng = v[5]
            }
        }
    };
};