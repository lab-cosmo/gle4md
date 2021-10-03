#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Load in pygtk and gtk

import gi
gi.require_version('Gtk', '3.0')
from gi. repository import Gtk
import gtk, gobject, cairo #, pango
import subprocess, math, re
from numpy import *

def logvalue(scale, value, fmt= "%4.2e"):
   return str( fmt % (10**value))


def adjwrap(value, lower, upper, step_increment, page_increment, page_size):
    # wrapper to Gtk.Adjustment to avoid depredation warning
    return Gtk.Adjustment(value=value, lower=lower, upper=upper,  
        step_increment=step_increment, page_increment=page_increment, page_size=page_size)

# Define the main window
class deltablock:
    def setw0(self,adj):
        self.omega0=adj.get_value()
        self.cback()
    def setdw(self,adj):
        self.domega=adj.get_value()
        self.cback()
    def setg(self,adj):
        self.gamma=adj.get_value()
        self.cback()
    def sett(self,adj):
        self.temp=adj.get_value()
        self.cback()

    def __init__(self, cback):
        self.cback=cback
        self.omega0=0.0
        self.domega=0.0
        self.gamma=0.0
        self.temp=0.0
        self.box = Gtk.Table(n_rows=4, n_columns=3, homogeneous=False)
        label=Gtk.Label(label="ω_0= "); label.show(); self.box.attach(label,0,1,0,1,Gtk.AttachOptions.SHRINK,Gtk.AttachOptions.SHRINK)
        label=Gtk.Label(label="Δω/ω_0= "); label.show(); self.box.attach(label,0,1,1,2,Gtk.AttachOptions.SHRINK,Gtk.AttachOptions.SHRINK)
        label=Gtk.Label(label="γ/(Δω ω_0)= ");   label.show(); self.box.attach(label,0,1,2,3,Gtk.AttachOptions.SHRINK,Gtk.AttachOptions.SHRINK)
        label=Gtk.Label(label="T= "); label.show(); self.box.attach(label,0,1,3,4,Gtk.AttachOptions.SHRINK,Gtk.AttachOptions.SHRINK)
#        label=Gtk.Label(label="10^"); label.show(); self.box.attach(label,1,2,0,1,Gtk.AttachOptions.SHRINK,Gtk.AttachOptions.SHRINK)
#        label=Gtk.Label(label="10^"); label.show(); self.box.attach(label,1,2,1,2,Gtk.AttachOptions.SHRINK,Gtk.AttachOptions.SHRINK);
#        label=Gtk.Label(label="10^"); label.show(); self.box.attach(label,1,2,2,3,Gtk.AttachOptions.SHRINK,Gtk.AttachOptions.SHRINK);
#        label=Gtk.Label(label="10^"); label.show(); self.box.attach(label,1,2,3,4,Gtk.AttachOptions.SHRINK,Gtk.AttachOptions.SHRINK)

        self.w0a = adjwrap(self.omega0, -5.0, 6.0, 0.01, 1.0, 1.0)
        self.w0a.connect("value_changed", self.setw0)
        self.w0s = Gtk.HScale(adjustment=self.w0a)
        self.w0s.connect("format-value", logvalue)
        self.w0s.set_digits(2);  self.w0s.set_value_pos(Gtk.PositionType.LEFT);
        self.w0s.show()
        self.box.attach(self.w0s,2,3,0,1,Gtk.AttachOptions.FILL|Gtk.AttachOptions.EXPAND,Gtk.AttachOptions.SHRINK)
        self.dwa = adjwrap(self.domega, -5.0, 6.0, 0.01, 1.0, 1.0)
        self.dwa.connect("value_changed", self.setdw)
        self.dws = Gtk.HScale(adjustment=self.dwa)
        self.dws.connect("format-value", logvalue)
        self.dws.set_digits(2);  self.dws.set_value_pos(Gtk.PositionType.LEFT)
        self.dws.show()
        self.box.attach(self.dws,2,3,1,2,Gtk.AttachOptions.FILL|Gtk.AttachOptions.EXPAND,Gtk.AttachOptions.SHRINK)
        self.ga = adjwrap(self.gamma, -5.0, 6.0, 0.01, 1.0, 1.0)
        self.ga.connect("value_changed", self.setg)
        self.gs = Gtk.HScale(adjustment=self.ga)
        self.gs.connect("format-value", logvalue)
        self.gs.set_digits(2);        self.gs.set_value_pos(Gtk.PositionType.LEFT)
        self.gs.show()
        self.box.attach(self.gs,2,3,2,3,Gtk.AttachOptions.FILL|Gtk.AttachOptions.EXPAND,Gtk.AttachOptions.SHRINK)
        self.ta = adjwrap(self.temp, -5.0, 6.0, 0.01, 1.0, 1.0)
        self.ta.connect("value_changed", self.sett)
        self.ts = Gtk.HScale(adjustment=self.ta)
        self.ts.connect("format-value", logvalue)
        self.ts.set_digits(2);        self.ts.set_value_pos(Gtk.PositionType.LEFT)
        self.ts.show()
        self.box.attach(self.ts,2,3,3,4,Gtk.AttachOptions.FILL|Gtk.AttachOptions.EXPAND,Gtk.AttachOptions.SHRINK)
        self.box.show()

class plotarea(Gtk.Frame):
    def __init__(self, css=None, border_width = 0):
        self.data=[]; self.datax=[]
        self.bhw=1.0
        self.expression=""
        super().__init__()
        
        self.surface = None     
        self.vexpand = True; self.hexpand=True   
        self.set_border_width(border_width)

        self.area = Gtk.DrawingArea()
        self.add(self.area)

        self.connect("draw", self.on_draw)
        self.area.connect('configure-event', self.on_configure)

    def init_surface(self, area):
        # Destroy previous buffer
        if self.surface is not None:
            self.surface.finish()
            self.surface = None

        # Create a new buffer
        self.surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, area.get_allocated_width(), area.get_allocated_height())

    def on_draw(self, area, context):
        if self.surface is not None:
            context.set_source_surface(self.surface, 0.0, 0.0)            
            context.paint()
        else:
            print('Invalid surface')


    def on_configure(self, area, event, data=None): 
        self.redraw()
        return False
        
                
    def redraw(self): 
        self.init_surface(self.area)
        context = cairo.Context(self.surface)
        self.draw(context)
        self.surface.flush()
        self.area.queue_draw()


    def draw(self, cr):
        width, height = self.area.get_allocated_width(), self.area.get_allocated_height()
        cr.scale(width/10, -height/10)        
        cr.translate(5, -5)

        # background grid
        cr.rectangle(-5,-5,10,10)
        cr.set_source_rgb(1, 1, 1)
        cr.fill()
        cr.set_line_width(0.01); cr.set_source_rgb(0,0,0)
        for i in range(-5,5):
            cr.move_to(i,-5); cr.line_to(i,5); cr.stroke()
            cr.move_to(-5,i); cr.line_to(5,i); cr.stroke()
        cr.set_line_width(0.025); cr.set_source_rgb(0,0,0)
        cr.move_to(0,-5); cr.line_to(0,5);
        cr.move_to(5,0); cr.line_to(-5,0); cr.stroke();
         
        cr.set_line_width(0.025)
        cr.set_source_rgb(0,0,0)
        if len(self.datax)>0: cr.move_to(self.datax[0],0)

        expr="math.log10("+self.expression+")";
        expr=re.sub(r'bhw',str(self.bhw),expr)
        for i in range(1,len(self.datax)):
            #cr.line_to(self.datax[i],math.log10(10**self.datax[i]*0.5*self.bhw/math.tanh(10**self.datax[i]*self.bhw*0.5)) )
            try:
                xexpr=re.sub(r'X',"(10**"+str(self.datax[i])+")",expr)
                yexpr=eval(xexpr)
                cr.line_to(self.datax[i],float(yexpr))
            except: pass


        cr.stroke()

        if len(self.datax)>0:
            cr.set_line_width(0.05);
            for j in range(0,len(self.data[0])):
                if(j%6==0):
                    cr.set_source_rgb(1,0,0)
                elif (j%6==1):
                    cr.set_source_rgb(0,0,1)
                elif (j%6==2):
                    cr.set_source_rgb(0,0.5,0)
                elif (j%6==3):
                    cr.set_source_rgb(0,0,0)
                elif (j%6==4):
                    cr.set_source_rgb(0.5,0.5,0)
                elif (j%6==5):
                    cr.set_source_rgb(0,0.5,0,5)

                cr.move_to(self.datax[0],self.data[0][j])
                for i in range(1,len(self.datax)):
                    cr.line_to(self.datax[i],self.data[i][j])
                cr.stroke()

class mainwin:
    def valuechanged(self):
        if(self.fauto): self.compute_gle()

    def selchanged(self):
        self.plot.datax=[]; self.plot.data=[];
        for el in self.gledata:
            self.plot.datax.append(math.log10(el[0]))
            dd=[]
            if self.glesel["K"]: dd.append(math.log10(el[4]));
            if self.glesel["H"]: dd.append(math.log10(el[5]));
            if self.glesel["kv"]: dd.append(math.log10(el[2]/el[0]));
            if self.glesel["kk"]: dd.append(math.log10(el[3]/el[0]));
            if self.glesel["kh"]: dd.append(math.log10(el[1]/el[0]));
            if self.glesel["cqq"]: dd.append(math.log10(el[6]));
            if self.glesel["cpp"]: dd.append(math.log10(el[7]));
            if self.glesel["spectrum"]: dd.append(math.log10(math.fabs(el[13])));
            self.plot.data.append(dd)
        self.plot.bhw=self.bhw
        self.plot.expression=self.expression.get_text()
        self.plot.redraw()


    def set_auto(self, toggle, user1):
        self.fauto=self.auto.get_active()

    def set_sel(self, toggle, user1):
        self.glesel[user1]=toggle.get_active()
        self.selchanged()

    def setcppw0(self,adj):
        self.cppw0=10**adj.get_value()
        self.valuechanged()

    def apars(self,adj, user1):
        if (user1=="gamma0"): self.gamma0=adj.get_value();
        elif (user1=="temp0"): self.temp0=adj.get_value();
        self.valuechanged()

    def setx(self,adj):
        self.bhw=10**adj.get_value()
        self.valuechanged()

    def ddraw(self,adj):
        while len(self.deltas)>int(adj.get_value()):
            self.dbox.remove(self.deltas[len(self.deltas)-1].box); self.deltas.pop();
        while len(self.deltas)<int(adj.get_value()):
            self.deltas.append(deltablock(self.valuechanged))
            self.dbox.pack_start(self.deltas[len(self.deltas)-1].box,False,False,0)
        self.valuechanged()

    def __init__(self):
        self.fauto=False
        self.deltas=[]
        self.gledata=[]
        self.cppw0=1.0
        self.glesel={"K" : True, "H" : False, "kv" : True, "kh" : False, "kk" : False, "cpp" : False, "cqq" : True, "spectrum": False}
        self.gamma0=-20; self.temp0=0; self.bhw=1.0
        # Window and framework
        self.window = Gtk.Window(type=Gtk.WindowType.TOPLEVEL)
        self.window.set_title("Build and test a GLE as a sum of deltas.")
        self.window.set_border_width(10)
        self.window.connect("destroy", self.destroy)

        # A Button, with an action
        # Add it to the geometry
        # show the button
        h0box = Gtk.HBox(homogeneous=False, spacing=3)
        vbox = Gtk.VBox(homogeneous=False, spacing=2)
        hbox = Gtk.HBox(homogeneous=False, spacing=3)
        label=Gtk.Label(label="Output prefix")
        hbox.pack_start(label,False,False,0)
        label.show()
        self.prefix=Gtk.Entry()
        self.prefix.set_text("deltafit")
        hbox.pack_start(self.prefix,False,False,0)
        self.prefix.show()
        self.compute = Gtk.Button(label="Compute!")
        self.compute.connect("clicked", self.click_compute, None)
        hbox.pack_start(self.compute,False,False,0)
        self.compute.show()
        self.auto= Gtk.CheckButton(label="Auto-compute")
        self.auto.connect("toggled", self.set_auto, None)
        self.auto.show()
        hbox.pack_end(self.auto, True, True, 0)
        hbox.show()
        vbox.pack_start(hbox,False,False,0)

        hbox = Gtk.HBox(homogeneous=False, spacing=3)
        label=Gtk.Label(label="Plot expression")
        hbox.pack_start(label,False,False,0)
        label.show()
        self.expression=Gtk.Entry()
        self.expression.set_text("X*0.5*bhw/tanh(X*bhw*0.5)")
        hbox.pack_end(self.expression, True, True, 0)
        self.expression.show()
        hbox.show()
        vbox.pack_start(hbox,False,False,0)

        hbox = Gtk.HBox(homogeneous=False, spacing=3)
        label=Gtk.Label(label="b h w=")
        hbox.pack_start(label,False,False,0)
        label.show()
        self.xadj=adjwrap(1.0, -1.0, 3, 0.01, 1.0, 1.0)
        self.xadj.connect("value_changed", self.setx)
        self.xslider=Gtk.HScale(adjustment=self.xadj)
        self.xslider.connect("format-value", logvalue)
        self.xslider.set_digits(1)
        self.xslider.set_value_pos(Gtk.PositionType.LEFT)
        self.xslider.set_size_request(400,20)
        hbox.pack_start(self.xslider,True,True,0)
        self.xslider.show()
        hbox.show()
        vbox.pack_start(hbox,False,False,0)

        hbox = Gtk.HBox(homogeneous=False, spacing=3)
        check=Gtk.CheckButton(label="Cpp(ω)"); check.connect("toggled",self.set_sel,"spectrum"); check.show(); hbox.pack_start(check,False,False,0);
        label=Gtk.Label(label="ω0=")
        hbox.pack_start(label,False,False,0)
        label.show()
        self.xadj=adjwrap(0.0, -5.0, 6.0, 0.01, 1.0, 1.0) 
        self.xadj.connect("value_changed", self.setcppw0)
        self.xslider=Gtk.HScale(adjustment=self.xadj)
        self.xslider.connect("format-value", logvalue)
        self.xslider.set_digits(1)
        self.xslider.set_value_pos(Gtk.PositionType.LEFT)
        self.xslider.set_size_request(400,20)
        hbox.pack_start(self.xslider,True,True,0)
        self.xslider.show()
        hbox.show()
        vbox.pack_start(hbox,False,False,0)

        hbox = Gtk.HBox(homogeneous=False, spacing=3)
        label=Gtk.Label(label="n_δ=")
        hbox.pack_start(label,False,False,0)
        label.show()
        self.dadj=adjwrap(1.0, 1.0, 16.0, 1.0, 1.0, 1.0)
        self.dadj.connect("value_changed", self.ddraw)
        self.dslider=Gtk.HScale(adjustment=self.dadj)
        self.dslider.set_digits(0)
        self.dslider.set_value_pos(Gtk.PositionType.LEFT)
        self.dslider.set_size_request(400,20)
        hbox.pack_start(self.dslider,True,True,0)
        self.dslider.show()
        hbox.show()
        vbox.pack_start(hbox,False,False,0)

        hbox = Gtk.HBox(homogeneous=False, spacing=3)
        check=Gtk.CheckButton(label="K(ω)"); check.set_active(True); check.connect("toggled",self.set_sel,"K"); check.show(); hbox.pack_start(check,False,False,0);
        check=Gtk.CheckButton(label="H(ω)"); check.connect("toggled",self.set_sel,"H"); check.show(); hbox.pack_start(check,False,False,0)
        check=Gtk.CheckButton(label="κ_V"); check.set_active(True); check.connect("toggled",self.set_sel,"kv"); check.show(); hbox.pack_start(check,False,False,0)
        check=Gtk.CheckButton(label="κ_K"); check.connect("toggled",self.set_sel,"kk"); check.show(); hbox.pack_start(check,False,False,0)
        check=Gtk.CheckButton(label="κ_H"); check.connect("toggled",self.set_sel,"kh"); check.show(); hbox.pack_start(check,False,False,0)
        check=Gtk.CheckButton(label="<K>"); check.connect("toggled",self.set_sel,"cpp"); check.show(); hbox.pack_start(check,False,False,0)
        check=Gtk.CheckButton(label="<V>"); check.set_active(True); check.connect("toggled",self.set_sel,"cqq"); check.show(); hbox.pack_start(check,False,False,0)
        hbox.show()
        vbox.pack_start(hbox,False,False,0)

        hbox = Gtk.HBox(homogeneous=False, spacing=3)
        label=Gtk.Label(label="γ_0=")
        hbox.pack_start(label,False,False,0)
        label.show()
        adj=adjwrap(self.gamma0, -20.0, 6.0, 0.01, 1.0, 1.0)
        adj.connect("value_changed", self.apars,"gamma0")
        slider=Gtk.HScale(adjustment=adj)
        slider.connect("format-value", logvalue)
        slider.set_digits(2)
        slider.set_value_pos(Gtk.PositionType.LEFT)
        slider.set_size_request(400,20); slider.show()
        hbox.pack_start(slider,True,True,0); hbox.show()
        vbox.pack_start(hbox,False,False,0)

        hbox = Gtk.HBox(homogeneous=False, spacing=3)
        label=Gtk.Label(label="T_0=")
        hbox.pack_start(label,False,False,0)
        label.show()
        adj=adjwrap(self.temp0, -5.0, 6.0, 0.01, 1.0, 1.0)
        adj.connect("value_changed", self.apars,"temp0")
        slider=Gtk.HScale(adjustment=adj)
        slider.connect("format-value", logvalue)
        slider.set_digits(2)
        slider.set_value_pos(Gtk.PositionType.LEFT)
        slider.set_size_request(400,20); slider.show()
        hbox.pack_start(slider,True,True,0); hbox.show()
        vbox.pack_start(hbox,False,False,0)


        self.dbox=Gtk.VBox(homogeneous=False, spacing=3)
        vbox.pack_start(self.dbox,False,False,0)
        self.dbox.show()
        hbox.show()
        self.plot=plotarea()
        
        vbox.show()
        h0box.pack_start(vbox,False,False,0)
        h0box.pack_end(self.plot, True, True, 0)
        h0box.show()
        self.window.add(h0box)
        self.plot.show()
        self.plot.set_size_request(500,500)

        # Show the window
        self.window.show_all()
        self.ddraw(self.dadj)
        self.compute_gle()

# Callback function for use when the button is pressed
    def click_compute(self, but, user1):
        self.compute_gle()

    def compute_gle(self):
        nd=int(len(self.deltas));ne=2*nd+1;
        amat=[[1e-20]*ne for i in range(ne)];
        dmat=[[1e-20]*ne for i in range(ne)];
        amat[0][0]=10**self.gamma0; dmat[0][0]=2*10**self.gamma0*10**self.temp0;
        for i in range(0,nd):
            w0=10**self.deltas[i].omega0; dw=10**self.deltas[i].domega*w0;
            g=math.sqrt(10**self.deltas[i].gamma*dw)

            #amat[0][2*i+1]=amat[0][2*i+2]=g; 
            #amat[2*i+1][0]=amat[2*i+2][0]=-g;             
            #amat[2*i+1][2*i+1]=amat[2*i+2][2*i+2]=dw
            #amat[2*i+1][2*i+2]=w0
            #amat[2*i+2][2*i+1]=-w0
            #dmat[2*i+1][2*i+1]=dmat[2*i+2][2*i+2]=2*dw*10**self.deltas[i].temp;
            amat[2*i+1][2*i+1]=dw;
            amat[2*i+1][2*i+2]=w0;
            amat[2*i+2][2*i+1]=-w0;
            amat[0][2*i+1]=g
            amat[2*i+1][0]=-g
            dmat[2*i+1][2*i+1]=2*dw*10**self.deltas[i].temp;

        afilename=self.prefix.get_text()+".a"
        dfilename=self.prefix.get_text()+".d"
        ofilename=self.prefix.get_text()+".dat"
        afile=open(afilename,'w')
        dfile=open(dfilename,'w')
        ofile=open(ofilename,'w')

        afile.write("rows "+str(ne)+"\n")
        afile.write("cols "+str(int(self.dadj.get_value())*2+1)+"\n")
        afile.write("data "+str( (int(self.dadj.get_value())*2+1)**2)+"\n")
        for i in range(0,ne):
            for j in range(0,ne):
                afile.write(str(amat[i][j])+" ")
            afile.write("\n")

        dfile.write("rows "+str(ne)+"\n")
        dfile.write("cols "+str(int(self.dadj.get_value())*2+1)+"\n")
        dfile.write("data "+str( (int(self.dadj.get_value())*2+1)**2)+"\n")
        for i in range(0,ne):
            for j in range(0,ne):
                dfile.write(str(dmat[i][j])+" ")
            dfile.write("\n")

        afile.close()
        dfile.close()        
        pgle=subprocess.Popen(["gle-analyze","-wi","1e-5","-wf","1e5","-np","500","-a",afilename,"-d",dfilename,"-w0",str(self.cppw0)],stdout=ofile)
        pgle.wait()
        ofile.close()
        ofile=open(ofilename,'r')
        ofile.readline(); ofile.readline(); # skips header
        self.gledata=[]
        while 1:
            line=ofile.readline()
            if not line:
                break
            ldata=[]
            for el in line.split():
                ldata.append(float(el))
            self.gledata.append(ldata)
        
        self.selchanged()


# Destroy method causes appliaction to exit
# when main window closed

    def destroy(self, widget, data=None):
        Gtk.main_quit()

# All PyGTK applicatons need a main method - event loop

    def main(self):
        Gtk.main()

if __name__ == "__main__":
    base = mainwin()
    base.main()
