
/*************************************************************************
*
*************************************************************************/
class RealtimeClock {
    
    /*********************************************************************
    *
    *********************************************************************/
    constructor() {
        
        this.time0 = new Date().getTime()
        
    } // end of constructor
        
    /*********************************************************************
    *
    *********************************************************************/
    getTime() {
        
        return (new Date().getTime() - this.time0)/1000.0;
        
    } // end of getTime method
    
    
} // end of RealtimeClock class


/*************************************************************************
*
*************************************************************************/
class AdjustableClock {
    
    /*********************************************************************
    *
    *********************************************************************/
    constructor() {
        
        this.clock = new RealtimeClock();
        this.time0 = this.clock.getTime();
        this.time1 = this.time0;
        this.rate = 1.0;
        
    } // end of constructor
    

    /*********************************************************************
    *
    *********************************************************************/
    setRate(rate) {
        
        if(rate<1.0) rate = 1.0;
        
        console.log('clock rate=',rate);
        
        const now = this.clock.getTime();
        this.time1 =  this.time1 + (now-this.time0)*this.rate;
        this.time0 = now;
        this.rate = rate;
        
        //console.log(this.time0, this.time1, this.rate);
    }
        
    /*********************************************************************
    *
    *********************************************************************/
    getTime() {
        
        //console.log(this.time0, this.time1, this.rate);
        
        const now = this.clock.getTime();
        return this.time1 + (now-this.time0)*this.rate;
        
    } // end of getTime method
    
    
} // end of RealtimeClock class


/*************************************************************************
*
*************************************************************************/
class CompositeAccelerator {
     
    /*********************************************************************
    *
    *********************************************************************/
    constructor() {
        
        this.accelerators = [];
        
    } // end of constructor
    
     
    /*********************************************************************
    *
    *********************************************************************/
    add(accelerator) {
        this.accelerators.push(accelerator);
    }
    
     
    /*********************************************************************
    *
    *********************************************************************/
    acceleration(state) { 
        
        let ax = 0.0;
        let ay = 0.0;
        let az = 0.0;
        
        for(const acc of this.accelerators) {
            
            const a = acc.acceleration(state);
            ax += a[0];
            ay += a[1];
            az += a[2];
            
        }
        
        return [ax, ay, az];
        
    } // end of acceleration method
    
    /*********************************************************************
    *
    *********************************************************************/
    handleCollision(state) {
        
        let max_v_impact = 0.0;
        for(const acc of this.accelerators) {
            if(acc.handleCollision) {
                const v_impact= acc.handleCollision(state);
                if(v_impact > max_v_impact) {
                    max_v_impact = v_impact;
                }
            }
        }
        
        return max_v_impact;
        
    } // end of handleCollision method
    
} // end of CompositeAccelerator class

/*************************************************************************
*
*************************************************************************/
class State {
    
    /*********************************************************************
    *
    *********************************************************************/
    constructor() {
    
        this.time = 0.0;
        this.x = [0.0, 0.0, 0.0];
        this.v = [0.0, 0.0, 0.0];
        this.a = [0.0, 0.0, 0.0];

    } // end of constructor
    
    /*********************************************************************
    *
    *********************************************************************/    
    set(time, x, y, z, vx, vy, vz, accelerator) {
        
        this.time = time;
        
        this.x[0] = x;
        this.x[1] = y;
        this.x[2] = z;
        
        this.v[0] = vx;
        this.v[1] = vy;
        this.v[2] = vz;   
        
        if(accelerator) this.a = accelerator.acceleration(this);
        else            this.a = 0.0;
        
    }
    
    /*********************************************************************
    *
    *********************************************************************/    
    copyFrom(state) {
        
        this.time = state.time;
        
        this.x[0] = state.x[0];
        this.x[1] = state.x[1];
        this.x[2] = state.x[2];
        
        this.v[0] = state.v[0];
        this.v[1] = state.v[1];
        this.v[2] = state.v[2];        
        
        this.a[0] = state.a[0];
        this.a[1] = state.a[1];
        this.a[2] = state.a[2];  
        
    } // end of copyFrom method
    
        
    /*********************************************************************
    *
    *********************************************************************/
    toString() {
        
        return this.time+" ("+this.x[0]+", "+this.x[1]+", "+this.x[2]+") ("+this.v[0]+", "+this.v[1]+", "+this.v[2]+")";
    }

    
    
} // end of State class

/*************************************************************************
*
*************************************************************************/
class RungeKutta {
    
    /*********************************************************************
    *
    *********************************************************************/
    constructor() {
        
        this.order = 6;
        this.indep_coef = [0.0, 1.0/5.0, 3.0/10.0, 3.0/5.0, 1.0, 7.0/8.0];
        
        this.depen_coef = [[0.0,   0.0, 0.0, 0.0, 0.0, 0.0],
                           [1.0/5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [3.0/40.0, 9.0/40.0, 0.0, 0.0, 0.0, 0.0],
                           [3.0/10.0, -9.0/10.0, 6.0/5.0, 0.0, 0.0, 0.0],
                           [-11.0/54.0, 5.0/2.0, -70.0/27.0, 35.0/27.0, 0.0, 0.0],
                           [1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0, 0.0]];
                           
        this.result_coef = [37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0];
                           
        this.error_coef = [this.result_coef[0] - 2825.0/27648.0,
                           this.result_coef[1] - 0.0,
                           this.result_coef[2] - 18575.0/48384.0,
                           this.result_coef[3] - 13525.0/55296.0,
                           this.result_coef[4] - 277.0/14336.0,
                           this.result_coef[5] - 1.0/4.0];
        
    } // end of constructor
    
    /*********************************************************************
    *
    *********************************************************************/
    tryAdvance(state, delta_time, accelerator) {
        
        const states = [state];
        for(let i=1; i<this.order; ++i) {
            states.push(new State());
        }
        
        /**************************
        * take intermediate steps *
        **************************/
        for(let i=1; i<this.order; ++i) {

            const new_time = state.time + delta_time*this.indep_coef[i];

            let new_x  = state.x[0];
            let new_y  = state.x[1];
            let new_z  = state.x[2];

            let new_vx = state.v[0];
            let new_vy = state.v[1];
            let new_vz = state.v[2];


            for(let j=0; j<i; ++j) {

                new_x  += states[j].v[0] * this.depen_coef[i][j] * delta_time;
                new_y  += states[j].v[1] * this.depen_coef[i][j] * delta_time;
                new_z  += states[j].v[2] * this.depen_coef[i][j] * delta_time;

                new_vx += states[j].a[0] * this.depen_coef[i][j] * delta_time;
                new_vy += states[j].a[1] * this.depen_coef[i][j] * delta_time;
                new_vz += states[j].a[2] * this.depen_coef[i][j] * delta_time;
            }

            states[i].set(new_time, 
                          new_x,  new_y,  new_z,
                          new_vx, new_vy, new_vz,
                          accelerator);

        } // end of outer loop

        /***********************
        * get the final result *
        ***********************/
        let new_x  = state.x[0];
        let new_y  = state.x[1];
        let new_z  = state.x[2];

        let new_vx = state.v[0];
        let new_vy = state.v[1];
        let new_vz = state.v[2];

        let error_x = 0.0;
        let error_y = 0.0;
        let error_z = 0.0;

        for(let j=0; j<this.order; ++j) {

            new_x += states[j].v[0] * this.result_coef[j] * delta_time;
            new_y += states[j].v[1] * this.result_coef[j] * delta_time;
            new_z += states[j].v[2] * this.result_coef[j] * delta_time;

            new_vx += states[j].a[0] * this.result_coef[j] * delta_time;
            new_vy += states[j].a[1] * this.result_coef[j] * delta_time;
            new_vz += states[j].a[2] * this.result_coef[j] * delta_time;

            error_x += states[j].v[0] * this.error_coef[j] * delta_time;
            error_y += states[j].v[1] * this.error_coef[j] * delta_time;
            error_z += states[j].v[2] * this.error_coef[j] * delta_time;

        }
        
        const result = new State();
        result.set(state.time + delta_time,
                   new_x, new_y, new_z,
                   new_vx, new_vy, new_vz,
                   accelerator);
        
        const error = Math.sqrt(error_x*error_x + 
                                error_y*error_y + 
                                error_z*error_z);
        
        return [result, error];

    } // end of tryAdvance method
    
    /*********************************************************************
    *
    *********************************************************************/    
    advance(state, new_time, accelerator, accuracy) {
        
        /*******************************
        * try the whole step in one go *
        *******************************/
        let time0 = state.time;
        let delta_time = new_time - time0;
        const [result, error] = this.tryAdvance(state, delta_time, accelerator);
        
        if(error>accuracy) {
            /************************
            * subdivide recursively *
            ************************/
            console.log('subdivide');
            const target = delta_time * Math.pow(accuracy/error, 0.2)*.9;
            let nsteps = Math.ceil(Math.pow(accuracy/error, 0.2)*0.9);
            if(nsteps<2) nsteps = 2;
            
            delta_time = delta_time/nsteps;
            for(let step=0; step<nsteps; ++step) {
                this.advance(state, time0 + delta_time*(step+1), 
                             accelerator, accuracy);
            }
            
        } else {
            /***********************
            * we made it in one go *
            ***********************/
            state.copyFrom(result);
            
        }
            
    } // end of advance method
    
    
} // end of RungeKutta class
