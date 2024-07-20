/*************************************************************************
*
*************************************************************************/
class StationaryPath {
    
    /*********************************************************************
    *
    *********************************************************************/ 
    constructor(x, y, z) { 
        this.x = [x, y, z];
    }
    
    getPosition(time) {
        return this.x;
    }
    
    getVelocity(time) {
        return [0.0, 0.0, 0.0];
    }
    
} // end of StationaryPath class

/*************************************************************************
*
*************************************************************************/
class CircularPath {
    
    /*********************************************************************
    *
    *********************************************************************/ 
    constructor(radius, period, phase, x, y, z) { 
        this.radius = radius;
        this.center = [x, y, z];
        this.speed=2*Math.PI/period;
        this.phase = phase*Math.PI/180.0;
    }
    
    getPosition(time) {
        
        const angle = this.speed*time + this.phase;
        
        const x = this.center[0] + this.radius*Math.cos(angle);
        const y = this.center[1] + this.radius*Math.sin(angle);
        const z = this.center[2];
        
        return [x, y, z];
    }
    
    getVelocity(time) {
        
        const angle = this.speed*time + this.phase;
        
        const vx = -this.radius*this.speed*Math.sin(angle);
        const vy =  this.radius*this.speed*Math.cos(angle);
        const vz = 0.0;
        
        return [vx, vy, vz];

    }
    
} // end of CircularPath class

/*************************************************************************
*
*************************************************************************/
class Body {
    
    /*********************************************************************
    *
    *********************************************************************/ 
    constructor(name, mass, radius, path) {
        
        this.name = name;
        this.mass = mass;
        this.radius = radius;
        this.path = path;
        
        this.G = 6.67430e-11;
        this.GM = this.G * this.mass;
        
    } // end of constructor
    
    /*********************************************************************
    *
    *********************************************************************/    
    acceleration(state) {
        
        const pos = this.path.getPosition(state.time);
        
        const dx = state.x[0] - pos[0];
        const dy = state.x[1] - pos[1];
        const dz = state.x[2] - pos[2];

        const r2 = dx*dx + dy*dy + dz*dz;
        const total = this.GM/r2;
        const r = Math.sqrt(r2);
        
        return [-dx/r*total, -dy/r*total, -dz/r*total];
        
    } // end of gravity method
    
    /*********************************************************************
    *
    *********************************************************************/     
    handleCollision(state) {
        
        const pos = this.path.getPosition(state.time);
        
        const dx = state.x[0] - pos[0];
        const dy = state.x[1] - pos[1];
        const dz = state.x[2] - pos[2];
        
        const dist2 = dx*dx + dy*dy + dz*dz
        
        const alt = Math.sqrt(dist2) - this.radius;
        
        state.v_impact = 0.0;
        
        if(alt >= 0.0) return 0.0;

        /*******************************************************
        * backtrack position along the relative velocity vector
        ********************************************************/
        const vel = this.path.getVelocity(state.time);
        const vx = state.v[0] - vel[0];
        const vy = state.v[1] - vel[1];
        const vz = state.v[2] - vel[2];
        
        
        const dp2 = dist2 - this.radius*this.radius;
        const v2 = vx*vx + vy*vy + vz*vz;
        const pv = dx*vx + dy*vy + dz*vz;
        
        const backtrack1 = (pv + Math.sqrt(pv*pv - dp2*v2))/v2;
        const backtrack2 = (pv - Math.sqrt(pv*pv - dp2*v2))/v2;
        
        let backtrack = backtrack1;
        if(Math.abs(backtrack2) < Math.abs(backtrack1)) {
            backtrack = backtrack2;
        }
        
        //console.log('backtrack=',backtrack)
        
        state.x[0] -= backtrack*vx;
        state.x[1] -= backtrack*vy;
        state.x[2] -= backtrack*vz;
        
        state.v[0] = vel[0];
        state.v[1] = vel[1];
        state.v[2] = vel[2];

        return Math.sqrt(v2);
            
    } // end of handleCollision method
    
    
} // end of Body class

/*************************************************************************
*
*************************************************************************/
class BodySystem {
    
    /*********************************************************************
    *
    *********************************************************************/
    constructor() {
        this.bodies = [];
    }
    
    /*********************************************************************
    *
    *********************************************************************/
    add(body) {
        this.bodies.push(body);
    }
    
    /*********************************************************************
    * Return the planet with the strongest gravity
    *********************************************************************/
    getDominant(state) {
        
        let max = 0.0;
        let dominant = null;
        for(const body of this.bodies) {
            const a = body.acceleration(state);
            const g2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
            
            if(g2 > max) {
                max = g2;
                dominant = body;
            }
        }
        
        return dominant;
            
        
    } // end of getDominant method
    
    
    /*********************************************************************
    *
    *********************************************************************/
    getBodyToTrack(time, width, height, centerx, centery) {
        
        for(const body of this.bodies) {
            
            const length = width + height + body.radius;
            const length2 = length*length;
            
            const pos = body.path.getPosition(time);
            const dx = pos[0] - centerx;
            const dy = pos[1] - centery;
            const distance2 = dx*dx + dy*dy;
            
            if(distance2 < length2) {
                return body;
            }
        } // end of loop over planets
        
        /********************************************
        * if we get here, nothing met our criteria
        * so return the default planet
        ********************************************/
        return this.bodies[0];
        
        
    } // end of getPlanetToTrack method
    
} // end of BodySystem class


/*************************************************************************
*
*************************************************************************/
class Atmosphere {
    

    /*********************************************************************
    *
    *********************************************************************/
    constructor(rocket, body, density, scale_height) {
        
        this.rocket = rocket;
        this.body = body;
        this.density0 = density;
        this.scale_height = scale_height;
        
    }
    

    /*********************************************************************
    *
    *********************************************************************/
    density(state) {
        
        if(this.density0 == 0.0) return 0.0;
        
        const pos_body = this.body.path.getPosition(state.time);
        
        const dx = state.x[0] - pos_body[0];
        const dy = state.x[1] - pos_body[1];
        const dz = state.x[2] - pos_body[2];
        

        
        let alt = Math.sqrt(dx*dx + dy*dy + dz*dz) - this.body.radius;
        if(alt<0) alt=0.0;
        //console.log('alt=',alt);
        
        return this.density0*Math.exp(-alt/this.scale_height);
        
    }
    
    
    /*********************************************************************
    *
    *********************************************************************/    
    acceleration(state) {
        
        const v0 = this.body.path.getVelocity(state.time);
        
        const vx = state.v[0] - v0[0];
        const vy = state.v[1] - v0[1];
        const vz = state.v[2] - v0[2];
        
        const vmag = Math.sqrt(vx*vx + vy*vy + vz*vz);
        
        //console.log('vmag=',vmag, 'vx=', vx,'vy=',vy);
        
        const density = this.density(state);
        
        const area= this.rocket.getStage().area;
        const mass = this.rocket.getMass(state.time);
        const drag = -vmag*this.density(state)*area/mass;
        
        return [drag*vx,
                drag*vy,
                drag*vz];
            
    }
            
} // end of Atmosphere class

/*************************************************************************
*
*************************************************************************/
class Rocket {
    
    /*********************************************************************
    *
    *********************************************************************/ 
    constructor(accelerator, x, y, vx, vy) { 

        this.state = new State();
        this.had_collision = true;
        
        this.rk = new RungeKutta();
        this.accelerator = new CompositeAccelerator();
        this.accelerator.add(accelerator);
        this.accelerator.add(this);
        
        this.stages = [];
        this.throttle = 0.0;
        this.headingx = 0.0;
        this.headingy = -1.0;
        
        const step = 1.0; // degrees
        this.sin_step = Math.sin(step/180.0*Math.PI);
        this.cos_step = Math.cos(step/180.0*Math.PI);
        
        /***********************
        * initialize the state *
        ***********************/
        this.state.set(0.0, x, y, 0.0, vx, vy, 0.0, this.accelerator);
        
    } // end of constructor
    
    
    
    /*********************************************************************
    *
    *********************************************************************/ 
    setCallback(callback) {
        this.callback = callback;
        
    }

    
    /*********************************************************************
    *
    *********************************************************************/    
    addStage(stage) {
        this.stages.push(stage);
        
    }
    
    /*********************************************************************
    *
    *********************************************************************/
    getStage() {
        
        return this.stages[this.stages.length -1];
        
    }
    
    /*********************************************************************
    *
    *********************************************************************/
    getStageCount() {
        
        return this.stages.length;
        
    }  
    
    /*********************************************************************
    *
    *********************************************************************/    
    getMass(time) {
        
        let mass = 0.0;
        for(const stage of this.stages) {
            mass += stage.payload + stage.fuel;
        }
        
        if(this.burn) {
            mass -= this.burn.burnt(time);
        }
        
        return mass;
        
    } // end of getMass method
    
    /*********************************************************************
    *
    *********************************************************************/ 
    setThrottle(throttle) { 
        
        this.throttle = throttle;
    }
    
    /*********************************************************************
    *
    *********************************************************************/    
    resetThrottle() {
        
        this.setThrottle(0.0);
        
        if(this.throttle_callback) {
            this.throttle_callback();
        }
    }
    
    /*********************************************************************
    *
    *********************************************************************/
    setThrottleResetCallback(callback) {
        
        this.throttle_callback = callback;
        
    }
    
    /*********************************************************************
    *
    *********************************************************************/ 
    rotateRight() {
        
        const x = this.headingx*this.cos_step - this.headingy*this.sin_step;
        const y = this.headingx*this.sin_step + this.headingy*this.cos_step;
        
        this.headingx = x;
        this.headingy = y;
        this.notify();
    }
    
    /*********************************************************************
    *
    *********************************************************************/ 
    rotateLeft() {
        
        const x =  this.headingx*this.cos_step + this.headingy*this.sin_step;
        const y = -this.headingx*this.sin_step + this.headingy*this.cos_step;
        
        this.headingx = x;
        this.headingy = y;
        this.notify();
    } 
    
    /*********************************************************************
    *
    *********************************************************************/ 
    thrust(time, rate) {
        const stage = this.getStage();
        const vexhaust = stage.vexhaust;
        
        return rate*vexhaust/this.getMass(time);
        
    }
    
    /*********************************************************************
    *
    *********************************************************************/ 
    acceleration(state) {
        
        if(!this.burn) {
            return [0.0, 0.0, 0.0];            
            
        } else {
            const vexhaust = this.getStage().vexhaust;
            const thrust   = this.thrust(state.time, this.burn.rate);

            return [this.headingx*thrust, 
                    this.headingy*thrust,
                    0.0];
        }
        
    } // end of acceleration method
    
    /*********************************************************************
    *
    *********************************************************************/ 
    notify() {
        
        if(this.callback) {
            this.callback(this);
        }
        
    } // end of notify
    
    /*********************************************************************
    *
    *********************************************************************/ 
    update(time) {
        
        //console.log('burn over', time-this.state.time, this.state.time, time);
        
        /*********************************
        * initialize the fuel accounting *
        *********************************/
        const stage = this.getStage();
        this.burn = stage.createBurn(this.state.time, time, this.throttle);
        
        /**********************
        * integrate the state *
        **********************/
        //const start = this.real_time.getTime();
        this.rk.advance(this.state, time, this.accelerator, 0.01);

        /********************
        * handle collisions *
        ********************/
        if(this.accelerator.handleCollision) {
            const v_impact = this.accelerator.handleCollision(this.state);
            if(v_impact > 0.0) {
                if(!this.had_collision) {
                    const mass = this.getMass(time);
                    if(mass*Math.abs(v_impact) >= this.getStage().crash_momentum) {
                        alert('Crash!');
                    }
                }
                    
                this.had_collision = true;
            } else {
                this.had_collision = false;
            }
            
        } // 
        
        /*********************
        * decrement the fuel *
        *********************/
        stage.burn(this.burn.burnt(time));
        this.burn = null;
        
        /******************************
        * stage if we are out of fuel *
        ******************************/
        if(this.getStage().fuel == 0.0 && this.stages.length > 1) {
            this.stages.pop();
            this.resetThrottle();
        }
            
        
    } // end of update method
    
    
} // end of Rocket class


/*************************************************************************
*
*************************************************************************/
class Stage {
    
    /*********************************************************************
    *
    *********************************************************************/ 
    constructor(payload, fuel, vexhaust, peak_flow, area, crash_momentum) { 
        
        this.payload = payload;
        this.max_fuel = fuel;
        this.fuel = fuel;
        this.vexhaust = vexhaust;
        this.peak_flow = peak_flow;
        this.area = area;
        this.crash_momentum = crash_momentum;

    }
    
    /*********************************************************************
    *
    *********************************************************************/
    createBurn(time0, time1, throttle) {
        
        let rate = this.peak_flow*throttle;
        if(rate*(time1-time0) > this.fuel) {
            rate = this.fuel/(time1-time0);
        }
        
        return new Burn(time0, rate);
        
        
    } // end of createBurn method
    
    
    /*********************************************************************
    *
    *********************************************************************/
    burn(fuel) {
        
        this.fuel -= fuel;
        
        if(this.fuel<1e-6) {
            this.fuel = 0.0;
        }
    }
    
} // end of stage class

/*************************************************************************
*
*************************************************************************/
class Burn {
       
    /*********************************************************************
    *
    *********************************************************************/ 
    constructor(start, rate) {
        
        this.start = start;
        this.rate = rate;
        
    } // end of constructor
    
       
    /*********************************************************************
    *
    *********************************************************************/ 
    burnt(time) {
        
        return this.rate*(time-this.start);
        
    }
        
} // end of Burn class

/*************************************************************************
*
*************************************************************************/
class Orbit {
    
    /*********************************************************************
    *
    *********************************************************************/ 
    constructor(state, body) {
        
        this.body = body;
        
        const time = state.time;
        
        const x_body = body.path.getPosition(time);
        const v_body = body.path.getVelocity(time);
        
        this.focus_x = x_body[0];
        this.focus_y = x_body[1];
        
        const x = state.x;
        const v = state.v;
        
        const vx = v[0] - v_body[0];
        const vy = v[1] - v_body[1];
        const vz = v[2] - v_body[2];
        
        const v2 = vx*vx + vy*vy + vz*vz;
        
        const dx = x[0] - x_body[0];
        const dy = x[1] - x_body[1];
        const dz = x[2] - x_body[2];
        const r = Math.sqrt(dx*dx + dy*dy + dz*dz);
        
        this.semi_major =  1.0/(2.0/r - v2/body.GM);
        
        /***************
        * eccentricity *
        ***************/
        const  lx = dy*vz - dz*vy;
        const  ly = dz*vx - dx*vz;
        const  lz = dx*vy - dy*vx;

        const  l2 = lx*lx + ly*ly + lz*lz;
        
        this.latus_rectum = l2/body.GM;
        this.semi_minor = Math.sqrt(this.latus_rectum * this.semi_major)

        this.eccentricity = Math.sqrt(1.0 - this.latus_rectum/this.semi_major);
        
        
        console.log('lr=',this.latus_rectum, 'e=',this.eccentricity);
        
        /***************************
        * orientation of the orbit *
        ***************************/
        const a = this.semi_major;
        const e = this.eccentricity;
        
        let delta = 0.0;
        
        if(this.latus_rectum < 1e-3) {
            delta =   -Math.PI;
        } else {
            delta = Math.acos((this.latus_rectum/r-1)/e);
        }

        // Now should this be positive or negative?
        const r_dot_v = dx*vx +  dy*vy +  dz*vz;
        const r_cross_v = dx*vy - dy*vx;
        if(r_dot_v * r_cross_v < 0. ) delta = - delta;

        // now find the absolute angle
        this.angle = 180.0/Math.PI*(Math.atan2(dy,dx) - delta);  
        
        console.log('delta=', delta, 'dx', dx, 'dy', dy);
        
        //console.log(this.semi_major, this.eccentricity, this.angle);
        
    } // end of constructor
    
    /*********************************************************************
    *
    *********************************************************************/ 
    angleToPoint(theta) {
        
        const cos = Math.cos(Math.PI/180.*(theta - this.angle));
        
        let r = 0.0
        if(this.latus_rectum < 1.0) {
            console.log('cos', cos, 'cos+1', cos+1, 'theta', theta, 'angle', this.angle);
            if(cos == -1.0) {
                r =  2.0*this.semi_major;
            } else {
                r =  0.0;
            }
            
        } else {

            r = this.latus_rectum/(1.0 + this.eccentricity*cos);

        }
        
        const rad = theta/180.0*Math.PI;

        const x = r*Math.cos(rad) + this.focus_x;
        const y = r*Math.sin(rad) + this.focus_y

        return [x, y];
        
        
    } // end of angleToPoint method
    
    
    /*********************************************************************
    *
    *********************************************************************/ 
    angleRange(r) {
        
        const costheta = (this.latus_rectum/r-1.0)/this.eccentricity;
        const theta = Math.abs(180.0/Math.PI*(Math.acos(costheta)));
        
        const angle1 = -theta + this.angle;
        const angle2 =  theta + this.angle;
        
        return [angle1, angle2];
        
    }
        
} // end of Orbit class
