/*************************************************************************
*
*************************************************************************/
class Plane {
    
    /*********************************************************************
    *
    *********************************************************************/
    constructor(width, height, system) {
        
//         this.width  = width;
//         this.height = height;
        this.system = system;
        
        /*******************************
        * initialize the body tracking *
        *******************************/
        this.tracked_body = system.bodies[0];     
        
        const pos = this.tracked_body.path.getPosition(0);
        this.centerx = pos[0]; // meters
        this.centery = pos[1]-this.tracked_body.radius; // meters
        this.scale = 10.0; // meters per pixel
        
        this.svg = document.createElementNS('http://www.w3.org/2000/svg','svg');
        
        
        
        this.svg.addEventListener('mousemove', (event) => {
            this.mouseMove(event.x, event.y, event.shiftKey || event.ctrlKey);
        });
        
        this.svg.addEventListener('mousedown', (event) => {
            this.mouseDown(event.x, event.y);
        });        
            
        this.svg.addEventListener('mouseup', (event) => {
            this.mouseUp(event.x, event.y);
        });  
        
        this.svg.addEventListener('dblclick', (event) => {
            if(this.center_callback) {
                this.center_callback();
            }
        }); 
        
        this.svg.setAttribute('tabindex', 0);

//         document.addEventListener('keypress', (event) => {
//             this.keyPressed(event.key);
//         });
        
        this.svg.setAttribute('width', '100%');
        this.svg.setAttribute('height', height);
        //this.svg.setAttribute('style', 'overflow: visible');


        this.items = []
        this.view_callbacks = [];
        
    } // end of constructor
    
    /*********************************************************************
    *
    *********************************************************************/ 
    onCenter(callback) {
        
        this.center_callback = callback;
    }
    
    
    /*********************************************************************
    *
    *********************************************************************/  
    addViewCallback(callback) {
        
        this.view_callback = callback; 
        
    }
    
    /*********************************************************************
    *
    *********************************************************************/
    notifyScaleChanged() {
        for(const callback of this.view_callbacks) {
            callback(this);
        }
    }
    
    /*********************************************************************
    *
    *********************************************************************/
    getWidth() {
        //console.log('this.svg.width=',this.svg.width.baseVal, this.svg.width.animVal);
        return this.svg.width.baseVal.value;

    }
    
    /*********************************************************************
    *
    *********************************************************************/
    getHeight() {
        return this.svg.height.baseVal.value;
        //console.log('height=',this.svg.getAttribute('height'));
    }
    
    /*********************************************************************
    *
    *********************************************************************/
    toPixels(x, y, time) {
        
        const pos = this.tracked_body.path.getPosition(time);
        
        const x1 = (x-this.centerx - pos[0])/this.scale + this.getWidth() /2;
        const y1 = (y-this.centery - pos[1])/this.scale + this.getHeight()/2;

        return [x1, y1];
        
    } // end of toPixels method
    
    /*********************************************************************
    *
    *********************************************************************/
    toMeters(x, y, time) {
        
        const pos = this.tracked_body.path.getPosition(time);
        
        const x1 = (x - this.getWidth() /2)*this.scale + this.centerx + pos[0];
        const y1 = (y - this.getHeight()/2)*this.scale + this.centery + pos[1];

        return [x1, y1];
        
    } // end of toMeters method 
    
    /*********************************************************************
    *
    *********************************************************************/
    distToPixels(r) {
        
        return r/this.scale;

    } // end of distToPixels method  
    
    
    /*********************************************************************
    *
    *********************************************************************/
    distToMeters(r) {
        
        //console.log('dist to meters r=',r,' scale=',this.scale);
        
        return r*this.scale;

    } // end of distToMeters method    
    
        
    /*********************************************************************
    *
    *********************************************************************/
    velocityToPixels(state) {
        
        const vb = this.tracked_body.path.getVelocity(state.time);
        
        const vx = state.v[0] - vb[0];
        const vy = state.v[1] - vb[1];
        const vz = state.v[2] - vb[2];
        
        const vv = Math.sqrt(vx*vx + vy*vy + vz*vz);
        
        console.log('state=',state);
        console.log('vb=',vb);
        console.log('delta v=',vx, vy, vz);
        console.log('vv=', vv);
        
        return this.distToPixels(vv);
        
    } // end of velocityToPixels method
   
    /*********************************************************************
    *
    *********************************************************************/
    setScale(scale) {
        
        this.scale = scale;
        this.notifyScaleChanged();
    }
    
    /*********************************************************************
    *
    *********************************************************************/ 
    zoomOut() {
        
        this.setScale(this.scale*1.5);
    }
    
    /*********************************************************************
    *
    *********************************************************************/ 
    zoomIn() {
        
        this.setScale(this.scale/1.5);
    }
    
    /*********************************************************************
    *
    *********************************************************************/
    centerOn(x, y, time) {
        
        const width = this.getWidth();
        const height = this.getHeight();
        const pos = this.toMeters(width/2, height/2, time);
        
//         console.log('width=',width);
//         console.log('height=',height);
//         console.log('center to meters', pos);
//         console.log('new center', x, y);
        
        this.centerx += x - pos[0];
        this.centery += y - pos[1];

    } // emd of centerOn method   
    
    /*********************************************************************
    *
    *********************************************************************/ 
    add(item) {
        
        item.setPlane(this);
        
        this.items.push(item);
        this.svg.appendChild(item.getElement());
        
    }
    
    /*********************************************************************
    *
    *********************************************************************/ 
    trackBody(time) {
        
        const width = this.distToMeters(this.getWidth());
        const height = this.distToMeters(this.getHeight());
        
        const pos = this.toMeters(this.getWidth(), this.getHeight(), time);
        
        
        const tracked_body = this.system.getBodyToTrack(time, width, height, 
                                                        pos[0], pos[1]);
    
        if(this.tracked_body != tracked_body) {
            console.log('tracking ', tracked_body.name);
            
            const old_pos = this.tracked_body.path.getPosition(time);
            const new_pos =      tracked_body.path.getPosition(time);
            
//             console.log('before', this.centerx, this.centery);
//             console.log('old_pos', old_pos);
//             console.log('new_pos', new_pos);
            this.centerx += old_pos[0];
            this.centerx -= new_pos[0];
            
            this.centery += old_pos[1];
            this.centery -= new_pos[1];
            
            //console.log('after', this.centerx, this.centery);
            
            this.tracked_body = tracked_body;

        }
        
    } // end of updateBodyTracking method
    
    
    /*********************************************************************
    *
    *********************************************************************/ 
    update(time) {
        
        this.trackBody(time);
        
        for(const item of this.items) {
            
            item.update(time);
        
        }
    } // end of update method
    
    /*********************************************************************
    *
    *********************************************************************/ 
    mouseDown(x, y) {
        
        this.dragging = true;
        this.mousex0 = x;
        this.mousey0 = y;
        
        this.centerx0 = this.centerx;
        this.centery0 = this.centery;

    }
    
    /*********************************************************************
    *
    *********************************************************************/ 
    mouseMove(x, y, stretch) {
        
        
        if(!this.dragging) return;


        const dx = x - this.mousex0;
        const dy = y - this.mousey0;
        
        this.centerx = this.centerx0 - this.distToMeters(dx);
        this.centery = this.centery0 - this.distToMeters(dy);

    } 
    
    
    /*********************************************************************
    *
    *********************************************************************/ 
    mouseUp(x, y) {
        
        this.dragging = false;
        
    }      
    
    
} // end of Plane class


/*************************************************************************
*
*************************************************************************/
class PlaneElement {
    
    /*********************************************************************
    *
    *********************************************************************/
    constructor(type) { 
        
        this.elem = document.createElementNS('http://www.w3.org/2000/svg',type);
    }
    
    /*********************************************************************
    *
    *********************************************************************/
    setPlane(plane) {
        this.plane = plane;
    }
    
    /*********************************************************************
    *
    *********************************************************************/
    getElement() { return this.elem; } 
    
    
    /*********************************************************************
    *
    *********************************************************************/
    update(time) {}    
    
} // end of PlaneElement super class


/*************************************************************************
*
*************************************************************************/
class GroupElement extends PlaneElement {
    
    /*********************************************************************
    *
    *********************************************************************/
    constructor() { 
        
        super('g');
        this.members = [];
    }
    
        
    /*********************************************************************
    *
    *********************************************************************/
    setPlane(plane) {
        
        super.setPlane(plane);
                
        for(const elem of this.members) {
            elem.setPlane(plane);
        }
        
    }
    
    /*********************************************************************
    *
    *********************************************************************/
    addElement(elem) {
        
        const group = this.getElement();
        group.appendChild(elem.getElement());
        
        this.members.push(elem);
        
        return elem;
        
    }
    
    
    /*********************************************************************
    *
    *********************************************************************/
    updateGroup(time) {
        
    }
    

    /*********************************************************************
    *
    *********************************************************************/
    update(time) {
        
        this.updateGroup(time);
        
        for(const elem of this.members) {
            elem.update(time);
        }
        
        
    } // end of update method  
    
} // end of PlaneElement super class


/*************************************************************************
*
*************************************************************************/
class LabelElement extends PlaneElement {
    
    /*********************************************************************
    *
    *********************************************************************/
    constructor() {
        
        super('text');
        
        this.dx = 0.0;
        this.dy = 0.0;
        this.text = ''
        
        const elem = this.getElement();
        elem.appendChild(document.createTextNode(this.text));
        
        elem.setAttribute('style', 'font: 12px sans-serif;');
    }
    
    /*********************************************************************
    *
    *********************************************************************/
    setText(string) { 
        
        this.text = string;
        
        const elem = this.getElement();
        while(elem.firstChild) {
            elem.removeChild(elem.firstChild);
        }
        elem.appendChild(document.createTextNode(this.text));
        
        
        
    }
    
    /*********************************************************************
    *
    *********************************************************************/
    setColor(color) {
        const elem = this.getElement();
        elem.setAttribute('fill', color);
    }
    
    /*********************************************************************
    * in meters
    *********************************************************************/
    setPosition(x, y) {
        
        this.x = x;
        this.y = y;
        
    }    
    
    /*********************************************************************
    *
    *********************************************************************/
    setOffset(dx, dy) {
        
        this.dx = dx;
        this.dy = dy;
        
    }
    
        
    /*********************************************************************
    *
    *********************************************************************/    
    updateContent() {
        
    }    
    
        
    /*********************************************************************
    *
    *********************************************************************/
    update(time) {
        
        this.updateContent();
   
        const xy = this.plane.toPixels(this.x, this.y, time);

        const elem = this.getElement();
        elem.setAttribute('x', xy[0] + this.dx);
        elem.setAttribute('y', xy[1] + this.dy);

    }
    
} // end of LabelElement class

/*************************************************************************
*
*************************************************************************/
class BodyElement extends PlaneElement {
    
    /*********************************************************************
    *
    *********************************************************************/
    constructor(body) { 
        
        super('path');
        
        this.body = body;
        
        const elem = this.getElement();
        elem.setAttribute('stroke', 'cyan');
        elem.setAttribute('fill-opacity', 1.0);
        elem.setAttribute('fill', 'cyan');
    }

    
    /*********************************************************************
    *
    *********************************************************************/ 
    update(time) {
        
        //console.log('moon draw time', time);

        const pos = this.body.path.getPosition(time);

        
        const elem = this.getElement();

        const angles = [];
        for(let i=0; i<360; ++i) {
            angles.push(i);
        }
        
//         let tmp = orbit.angle;
//         while(tmp<0.0) tmp += 360.0;
//         while(tmp>360.0) tmp -= 360.0;
//         angles.push(tmp);
//         
//         tmp += 180;
//         if(tmp > 360.0) tmp -= 360.0;
//         angles.push(tmp);
        
//        angles.sort((a, b) => { return a - b; });

        const r = this.plane.distToPixels(this.body.radius);

        let path = '';
        for(const theta of angles) {

            const radians = theta/180.0 * Math.PI;
            const x = pos[0] + this.body.radius * Math.cos(radians);
            const y = pos[1] + this.body.radius * Math.sin(radians);
            
            const p = this.plane.toPixels(x, y, time);
            
            
            if(path.length == 0) {
                path += 'M '+p[0]+' '+p[1];
            } else {
                path += ' A '+r+' '+r+' 0 0 1 '+
                        p[0]+' '+p[1];
            }
            
        }
        
        path += ' Z';

        
        elem.setAttribute('d', path);        
        
        

    }
    
} // end of BodyElement class


/*************************************************************************
*
*************************************************************************/
class RocketElement extends PlaneElement {
    
    /*********************************************************************
    *
    *********************************************************************/
    constructor(rocket) { 
        
        super('line');
        
        this.rocket = rocket;
        
        const elem = this.getElement();
        elem.setAttribute('stroke', 'black');
        elem.setAttribute('stroke-width', 3);
    
        
    
    } // end of constructor

    
    /*********************************************************************
    *
    *********************************************************************/
    update(time) {
        
        this.rocket.update(time);
        
        const pos = this.rocket.state.x;
        const xy = this.plane.toPixels(pos[0], pos[1], time);
        
        const elem = this.getElement();
        elem.setAttribute('x1', xy[0]);
        elem.setAttribute('y1', xy[1]);
        
        const length = 20*this.rocket.getStageCount();
        elem.setAttribute('x2', xy[0] + this.rocket.headingx*length);
        elem.setAttribute('y2', xy[1] + this.rocket.headingy*length);
        
        
    }
    
    
} // end of RocketElement class


/*************************************************************************
*
*************************************************************************/
class OrbitPathElement extends PlaneElement {
    
    /*********************************************************************
    *
    *********************************************************************/
    constructor() { 
        
        super('path');

        const elem = this.getElement();
        elem.setAttribute('stroke', 'red');
        elem.setAttribute('fill-opacity', 0.0);
        
    } // end of constructor
    
    
    
    /*********************************************************************
    *
    *********************************************************************/
    setOrbit(orbit) { this.orbit = orbit; }    
    
    
    /*********************************************************************
    *
    *********************************************************************/
    update(time) {
        
        const elem = this.getElement();
        //elem.setAttribute('stroke-width', 1.0);
        
        const orbit = this.orbit; 

        if(orbit.latus_rectum == 0) {
            /****************
            * straight line *
            ****************/
            let path = '';
            
            const p1m = orbit.angleToPoint(orbit.angle);
            const p2m = orbit.angleToPoint(orbit.angle+180.0);
            
            const p1 = this.plane.toPixels(p1m[0], p1m[1], time);
            const p2 = this.plane.toPixels(p2m[0], p2m[1], time);
        
            path += 'M '+p1[0]+' '+p1[1];
            path += ' L '+p2[0]+' '+p2[1];
            elem.setAttribute('d', path);
            
        } else if(orbit.eccentricity < 1.0) {
            /**********
            * ellipse *
            **********/
            const angles = [];
            for(let i=0; i<360; ++i) {
                angles.push(i);
            }
            
            let tmp = orbit.angle;
            while(tmp<0.0) tmp += 360.0;
            while(tmp>360.0) tmp -= 360.0;
            angles.push(tmp);
            
            tmp += 180;
            if(tmp > 360.0) tmp -= 360.0;
            angles.push(tmp);
            
            angles.sort((a, b) => { return a - b; });

            const semi_major = this.plane.distToPixels(orbit.semi_major);
            const semi_minor = this.plane.distToPixels(orbit.semi_minor);

            let path = '';
            for(const theta of angles) {

                const pm = orbit.angleToPoint(theta);
                const p = this.plane.toPixels(pm[0], pm[1], time);
                
                
                if(path.length == 0) {
                    path += 'M '+p[0]+' '+p[1];
                } else {
                    path += ' A '+semi_major+' '+semi_minor+' '+orbit.angle+' 0 1 '+
                            p[0]+' '+p[1];
                }
                
            }
            
            path += ' Z';

            
            elem.setAttribute('d', path);
            
        } else {
            /************
            * hyperbola *
            ************/
            const width  = this.plane.getWidth();
            const height = this.plane.getHeight();
            const center = this.plane.toMeters(width/2, height/2, time);
            
            
            
            const dx = Math.abs(center[0] - this.orbit.focus_x);
            const dy = Math.abs(center[1] - this.orbit.focus_y);
            
            const max_r = dx + dy + this.plane.distToMeters(width) +
                                    this.plane.distToMeters(height);
            
            const range = this.orbit.angleRange(max_r);
            const angle1 = range[0];
            const angle2 = range[1];
            
            //console.log('max_r=',max_r,'angle1=',angle1,'angle2=',angle2);

            const npoints = 101;
            const delta = (angle2-angle1)/(npoints-1);
            let path = '';
            for(let i=0; i<npoints; ++i) {
                
                const theta = angle1 + i*delta;
                const pm = orbit.angleToPoint(theta);
                const p = this.plane.toPixels(pm[0], pm[1], time);
                
                if(path.length == 0) {
                    path += 'M '+p[0]+' '+p[1];
                } else {
                    path += ' L '+p[0]+' '+p[1];
                }
                
            } // end of loop over points
            
            elem.setAttribute('d', path);
            
        } // end if it's a hyperbola
        
    } // end of update method
    
} // end of OrbitPathElement class

/*************************************************************************
*
*************************************************************************/
class OrbitDecoration extends GroupElement {
    
    /*********************************************************************
    *
    *********************************************************************/
    constructor(angle) { 
        super();
        
        this.angle = angle;
        
        this.tic_length = 10.0; //pixels;
        
        this.label = this.addElement(new LabelElement());
        this.label.setOffset(0, 12);
        
        this.tic   = this.addElement(new PlaneElement('line'));
        this.tic.getElement().setAttribute('stroke', 'red');
    }
    
    
    /*********************************************************************
    *
    *********************************************************************/
    setAngle(angle) {
        this.angle = angle;
        
    }
    
    /*********************************************************************
    *
    *********************************************************************/
    setOrbit(orbit) {
        
        this.orbit = orbit;
    }
    
    /*********************************************************************
    *
    *********************************************************************/
    updateGroup(time) {
        
        const dir = this.orbit.angle + this.angle;
        const opos_m = this.orbit.angleToPoint(dir);
        const opos = this.plane.toPixels(opos_m[0], opos_m[1], time);
        
        
        const bpos_m = this.orbit.body.path.getPosition(time);
        const bpos = this.plane.toPixels(bpos_m[0], bpos_m[1], time);
        
        const dx = opos[0] - bpos[0];
        const dy = opos[1] - bpos[1];
        const dz = 0.0;

        const alt = this.plane.distToMeters(Math.sqrt(dx*dx + dy*dy + dz*dz)) - this.orbit.body.radius;
        

        if(Math.abs(alt)<1000.0) this.label.setText(alt.toFixed(1)+'m');
        else                     this.label.setText(Number((alt/1000).toFixed(1)).toLocaleString()+'km');
        
        /***********
        * tic mark *
        ***********/

        const x1 = opos[0];
        const y1 = opos[1];
        
        const norm = Math.sqrt(dx*dx + dy*dy);
        if(norm == 0.0) {
            // we should handle this in a nicer way
            return;
        }
        const xhat = dx/norm;
        const yhat = dy/norm;
        
        const x2 = x1 - this.tic_length*xhat;
        const y2 = y1 - this.tic_length*yhat;

        const line = this.tic.getElement();
        //line.setAttribute('stroke-width', 1.0);
        
        line.setAttribute('x1', x1); 
        line.setAttribute('y1', y1);  
        line.setAttribute('x2', x2);  
        line.setAttribute('y2', y2);
        
        const label_pos = this.plane.toMeters(x2, y2, time);
        
        this.label.setPosition(label_pos[0], label_pos[1]);
        
    }
        
} // end of OrbitDecoration class


/*************************************************************************
*
*************************************************************************/
class OrbitElement extends GroupElement {
    
        
    /*********************************************************************
    *
    *********************************************************************/
    constructor(rocket, system) { 
        
        super();
        
        this.rocket = rocket;
        this.system = system;
        
        this.path    = this.addElement(new OrbitPathElement());
        this.perigee = this.addElement(new OrbitDecoration(  0.0));
        this.apogee  = this.addElement(new OrbitDecoration(180.0));
    }
        
    /*********************************************************************
    *
    *********************************************************************/
    updateGroup(time) {
        
        const body = this.system.getDominant(this.rocket.state);
        const orbit = new Orbit(this.rocket.state, body);
        
        this.path.setOrbit(orbit);
        this.perigee.setOrbit(orbit);
        
        /******************************************************************
        * a cheat for hyperbolas. Instead of removing the apogee marker, 
        * we turn it into a second perigee marker 
        ******************************************************************/
        if(orbit.eccentricity>=1) {
            this.apogee.setAngle(0.0);
        } else {
            this.apogee.setAngle(180.0);
        }
        
        this.apogee.setOrbit(orbit);

    } // end of updateGroup method

} // end of OrbitElement class


/*************************************************************************
*
*************************************************************************/
class ThrustLabel extends LabelElement {

    /*********************************************************************
    *
    *********************************************************************/
    constructor(rocket) {
        
        super();
        this.rocket = rocket

        this.setOffset(6, 0);
        
    } // end of constructor
    

    /*********************************************************************
    *
    *********************************************************************/
    updateContent() {
        
        const rocket = this.rocket;
        
        const x = rocket.state.x[0];
        const y = rocket.state.x[1];
        this.setPosition(x, y);
        
        const stage = rocket.getStage();
        const rate = rocket.throttle* stage.peak_flow;
        const thrust = this.rocket.thrust(rocket.state.time, rate)/9.81;
        this.setText('Thrust: '+thrust.toFixed(2));
        
    }
    
} // end of ThrustLabel class

/*************************************************************************
*
*************************************************************************/
class DragLabel extends LabelElement {

    /*********************************************************************
    *
    *********************************************************************/
    constructor(rocket, atmosphere) {
        
        super();
        this.rocket = rocket;
        this.atmosphere = atmosphere;

        this.setOffset(6, -14);
        
    } // end of constructor
    

    /*********************************************************************
    *
    *********************************************************************/
    updateContent() {
        
        const rocket = this.rocket;
        
        
        const x = rocket.state.x[0];
        const y = rocket.state.x[1];
        this.setPosition(x, y);
        
        const a = this.atmosphere.acceleration(rocket.state);
        const drag = Math.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])/9.81;
        
        if(drag >= 0.01) this.setText('Drag: '+drag.toFixed(2));
        else             this.setText('');
        
    }
    
} // end of DragLabel class


/*************************************************************************
*
*************************************************************************/
class VelocityLabel extends LabelElement {

    /*********************************************************************
    *
    *********************************************************************/
    constructor(rocket, system) {
        
        super();
        this.rocket = rocket;
        this.system = system;

        this.setOffset(6, -28);
        
    } // end of constructor
    

    /*********************************************************************
    *
    *********************************************************************/
    updateContent() {
        
        const rocket = this.rocket;
        const body   = this.system.getDominant(rocket.state);
        
        
        const x = rocket.state.x[0];
        const y = rocket.state.x[1];
        this.setPosition(x, y);
        
        const pb = body.path.getPosition(rocket.state.time);
        const px = rocket.state.x[0] - pb[0];
        const py = rocket.state.x[1] - pb[1];
        const pz = rocket.state.x[2] - pb[2];
        
        const r = Math.sqrt(px*px + py*py + pz*pz);
        const rx = px/r;
        const ry = py/r;
        const rz = pz/r;
        
        const alt = r - body.radius;
        
        const vr = rocket.state.v;
        const vb = body.path.getVelocity(rocket.state.time);
        
        const vx = vr[0] - vb[0];
        const vy = vr[1] - vb[1];
        const vz = vr[2] - vb[2];
        
        const v_radial = vx*rx + vy*ry + vz*rz;
        
        const v = Math.sqrt(vx*vx + vy*vy + vz*vz);
        
        const m = this.rocket.getMass(this.rocket.state.time);
        
        const momentum = m*v;

        if(alt > 1e4) {
            this.setText("");
        } else {
            if(v_radial >0.0) {
                this.setColor('black');
            } else {
                if(momentum < this.rocket.getStage().crash_momentum) {
                    this.setColor('green');
                } else {
                    this.setColor('red');
                }
            }
        
            this.setText('v: '+(v/1000.0*3600.0).toFixed(1)+' km/h');
        }

        
    } // end of update method
    
} // end of VelocityLabel class

/*************************************************************************
*
*************************************************************************/
class FuelIndicator extends PlaneElement {
    
    /*********************************************************************
    *
    *********************************************************************/
    constructor(rocket) { 
        
        super('line');
        
        this.rocket = rocket;
        
        const elem = this.getElement();
        elem.setAttribute('stroke', 'green');
        elem.setAttribute('stroke-width', 3);
    
        this.length = 100;
        this.offsetx = 10;
        this.offsety = 6;
    
    } // end of constructor

    
    /*********************************************************************
    *
    *********************************************************************/
    update(time) {
        
        
        const elem = this.getElement();
        
        
        const pos = this.rocket.state.x;
        const xy = this.plane.toPixels(pos[0], pos[1], time);
        
        const x = xy[0]+this.offsetx;
        const y = xy[1]+this.offsety;
        
        elem.setAttribute('x1', x);
        elem.setAttribute('y1', y);
        
        const stage = this.rocket.getStage();
        const length = this.length * stage.fuel/stage.max_fuel;

        elem.setAttribute('x2', x + length);
        elem.setAttribute('y2', y);
        
        
    } // end of update method
    
    
} // end of FuelIndicator class



/*************************************************************************
*
*************************************************************************/
class ClockDisplay {
    
    /*********************************************************************
    *
    *********************************************************************/
    constructor(clock) { 
        
        this.clock = clock;
        
        this.node = document.createTextNode('0d 00:00:00');
        
        setInterval(() => {this.update();}, 1000);
        
    } // end of constructor
    
    /*********************************************************************
    *
    *********************************************************************/
    getElement() { return this.node; }    
     
    /*********************************************************************
    *
    *********************************************************************/
    update() {
        
        const now = Math.round(this.clock.getTime());
        const seconds = now%60;
        const minutes = Math.floor(now/60)%60;
        const hours   = Math.floor(now/3600)%24;
        const days    = Math.floor(now/(3600*24));
        
        this.node.nodeValue = days+'d '+String(hours).padStart(2,'0')+':'+
                                        String(minutes).padStart(2,'0')+':'+
                                        String(seconds).padStart(2,'0')+' '
        
        
        
    } // end of update method
    
} // end of ClockDisplay class
    























































