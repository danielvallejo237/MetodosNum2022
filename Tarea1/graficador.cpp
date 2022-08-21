/*
Código para graficar funciones
Recibe una función mediante consola y el intervalo [a,b] correspondiente al
dominio de la función

Código escrito por Daniel Vallejo Aldana
*/

#include <bits/stdc++.h>
#include <cmath>
#include <cairo.h>
#include <gtk/gtk.h>
#include "./fparser/fparser.hh"
#include <cstring>
#include <fstream>

#define INF 100000000 //Lo que nosotros conocemos como nuestra cota de infinito
#define MAX_PIXEL_WIDTH 450
#define MAX_PIXEL_HEIGHT 450
#define MAX_POINTS 1000
#define EXPANSION 20
using namespace std;

class Data
{
public:
  double mx,mi,fmax,fmin;
  Data(double mx, double mi,double fmax, double fmin)
  {
    this->mx=mx;
    this->mi=mi;
    this->fmax=fmax;
    this->fmin=fmin;
  }
  Data()
  {
    this->mx=0;
    this->mi=0;
    this->fmax=0;
    this->fmin=0;
  }
  void add_values(double mx, double mi,double fmax, double fmin)
  {
    this->mx=mx;
    this->mi=mi;
    this->fmax=fmax;
    this->fmin=fmin;
  }
  ~Data(){} //Nuevamente el constructor vacío
};

class Info
{
    public:
        double f;
        double x;
        int index;
        Info(double f,double x, double index)
        {
            this->f=f;
            this->x=x;
            this->index=index;
        }
        ~Info(){}
        void print()
        {
            cout<<x<<" "<<f<<" "<<index<<endl;
        }
};
namespace utils
{
    vector<double> linspace(double a, double b, int N)
    {
        vector<double> P(N);
        double stepsize=(b-a)/(double)N;
        for(int i=0;i<N;i++) P[i]=a+(double)(i+1)*stepsize;
        return P; //Vector de puntos que contienen las evaluaciones de las funciones f(x)
    }
    vector<double> evaluate(vector<double> xaxis, FunctionParser fp, pair<Info,Info> *mxmi)
    {
        //Evaluar los vectores en los puntos
        int sz=xaxis.size();
        vector<double> Y;
        double vals[1];
        for(vector<double>::iterator it=xaxis.begin();it!=xaxis.end();++it)
        {
            vals[0]=*it;
            Y.push_back(fp.Eval(vals));
        }
        int i=0;
        for(vector<double>::iterator it=Y.begin(),it2=xaxis.begin();it!=Y.end(),it2!=xaxis.end();++it,++it2)
        {
            if(*it<mxmi->first.f)
            {
                mxmi->first.f=*it;
                mxmi->first.x=*it2;
                mxmi->first.index=i;
            }
            if(*it>mxmi->second.f)
            {
                mxmi->second.f=*it;
                mxmi->second.x=*it2;
                mxmi->second.index=i;
            }
            i++;
        }
        return Y;
    }
    vector<pair<int,int>> ConvertCoordinates(vector<double> X,vector<double> F, pair<Info,Info> *mxmi, double a, double b)
    {
        double fmax=mxmi->second.f;
        double fmin=mxmi->first.f;
        double xmax=b;
        double xmin=a;
        vector<pair<int,int>> nonzero;
        for (vector<double>::iterator it=X.begin()+1,it2=F.begin()+1;it!=X.end();++it,++it2)
        {
            int x=(int)round(((*it-a)/(b-a))*MAX_PIXEL_WIDTH);
            int y=(int)round(((fmax-*it2)/(fmax-fmin))*MAX_PIXEL_HEIGHT);
            nonzero.push_back(make_pair(x,y));
        }
    return nonzero;
    }
}
vector<pair<int,int>>  NZ;
Data Datos;
static void do_drawing(cairo_t *);

static gboolean on_draw_event(GtkWidget *widget, cairo_t *cr)
{
  do_drawing(cr);
  return FALSE;
}

static void do_drawing(cairo_t *cr)
{
  cairo_set_source_rgb(cr, 0.886275, 0.278431, 0.101961); //Un naranja ubuntu para graficar las cosas
  cairo_set_line_width(cr, 3.5);
  for(vector<pair<int,int>>::iterator it=NZ.begin();it!=NZ.end()-1;++it)
  {
      cairo_move_to(cr, (*it).first+EXPANSION, (*it).second+EXPANSION);
      cairo_line_to(cr, (*(it+1)).first+EXPANSION, (*(it+1)).second+EXPANSION);
  }
  cairo_stroke(cr);
  cairo_set_source_rgb(cr, 0, 0, 0); //La linea de los ejes es de color negro para resaltar
  cairo_set_line_width(cr, 1.5);
  //Linea sobre el eje x
  cairo_stroke(cr);
  cairo_set_font_size(cr,9);
  cairo_move_to(cr,0+EXPANSION,MAX_PIXEL_WIDTH+EXPANSION);
  cairo_line_to(cr,MAX_PIXEL_HEIGHT+EXPANSION,MAX_PIXEL_WIDTH+EXPANSION);
  cairo_move_to(cr,MAX_PIXEL_HEIGHT+EXPANSION,MAX_PIXEL_WIDTH+3*EXPANSION/4);
  char* char_arr;
  double num=Datos.mx;
  string st=to_string(num);
  char_arr = &st[0];
  cairo_show_text(cr,char_arr);
  cairo_move_to(cr,0+EXPANSION,MAX_PIXEL_WIDTH+3*EXPANSION/4);
  num=Datos.mi;
  string st1=to_string(num);
  char_arr=&st1[0];
  cairo_show_text(cr,char_arr);
  cairo_move_to(cr,MAX_PIXEL_HEIGHT/2+EXPANSION,MAX_PIXEL_WIDTH+3*EXPANSION/4);
  num=(Datos.mx+Datos.mi)/2;
  string st2=to_string(num);
  char_arr = &st2[0];
  cairo_show_text(cr,char_arr);
  cairo_move_to(cr,MAX_PIXEL_HEIGHT/4+EXPANSION,MAX_PIXEL_WIDTH+3*EXPANSION/4);
  num=Datos.mi+(Datos.mx-Datos.mi)/4;
  string st3=to_string(num);
  char_arr = &st3[0];
  cairo_show_text(cr,char_arr);
  cairo_move_to(cr,3*MAX_PIXEL_HEIGHT/4+EXPANSION,MAX_PIXEL_WIDTH+3*EXPANSION/4);
  num=Datos.mi+3*(Datos.mx-Datos.mi)/4;
  st=to_string(num);
  char_arr=&st[0];
  cairo_show_text(cr,char_arr);
  //Linea sobre el eje Y
  cairo_move_to(cr,0+EXPANSION,0+EXPANSION);
  cairo_line_to(cr,0+EXPANSION,MAX_PIXEL_HEIGHT+EXPANSION);
  cairo_move_to(cr,0+EXPANSION+3,0+EXPANSION);
  num=Datos.fmax;
  st=to_string(num);
  char_arr=&st[0];
  cairo_show_text(cr,char_arr);
  cairo_move_to(cr,0+EXPANSION+3,MAX_PIXEL_HEIGHT+EXPANSION);
  num=Datos.fmin;
  st=to_string(num);
  char_arr=&st[0];
  cairo_show_text(cr,char_arr);
  cairo_move_to(cr,0+EXPANSION+3,MAX_PIXEL_HEIGHT/2+EXPANSION);
  num=(Datos.fmax+Datos.fmin)/2;
  st=to_string(num);
  char_arr=&st[0];
  cairo_show_text(cr,char_arr);
  cairo_move_to(cr,0+EXPANSION+3,MAX_PIXEL_HEIGHT/4+EXPANSION);
  num=Datos.fmin+3*(Datos.fmax+Datos.fmin)/4;
  st=to_string(num);
  char_arr=&st[0];
  cairo_show_text(cr,char_arr);
  cairo_move_to(cr,0+EXPANSION+3,3*MAX_PIXEL_HEIGHT/4+EXPANSION);
  num=Datos.fmin+(Datos.fmax+Datos.fmin)/4;
  st=to_string(num);
  char_arr=&st[0];
  cairo_show_text(cr,char_arr);
  cairo_stroke(cr);

  //Creación de un pequeño grid para mejor visualización
  //Lineas en Horizontal
  cairo_set_line_width(cr, 0.5);
  cairo_move_to(cr,0+EXPANSION+1,0+EXPANSION);
  cairo_line_to(cr,MAX_PIXEL_WIDTH+EXPANSION-1,0+EXPANSION);
  cairo_move_to(cr,0+EXPANSION+1,MAX_PIXEL_HEIGHT+EXPANSION);
  cairo_line_to(cr,MAX_PIXEL_WIDTH+EXPANSION-1,MAX_PIXEL_HEIGHT+EXPANSION);
  cairo_move_to(cr,0+EXPANSION+1,MAX_PIXEL_HEIGHT/2+EXPANSION);
  cairo_line_to(cr,MAX_PIXEL_WIDTH+EXPANSION-1,MAX_PIXEL_HEIGHT/2+EXPANSION);
  cairo_move_to(cr,0+EXPANSION+1,MAX_PIXEL_HEIGHT/4+EXPANSION);
  cairo_line_to(cr,MAX_PIXEL_WIDTH+EXPANSION-1,MAX_PIXEL_HEIGHT/4+EXPANSION);
  cairo_move_to(cr,0+EXPANSION+1,3*MAX_PIXEL_HEIGHT/4+EXPANSION);
  cairo_line_to(cr,MAX_PIXEL_WIDTH+EXPANSION-1,3*MAX_PIXEL_HEIGHT/4+EXPANSION);
  cairo_move_to(cr,0+EXPANSION+1,7*MAX_PIXEL_HEIGHT/8+EXPANSION);
  cairo_line_to(cr,MAX_PIXEL_WIDTH+EXPANSION-1,7*MAX_PIXEL_HEIGHT/8+EXPANSION);
  cairo_move_to(cr,0+EXPANSION+1,3*MAX_PIXEL_HEIGHT/8+EXPANSION);
  cairo_line_to(cr,MAX_PIXEL_WIDTH+EXPANSION-1,3*MAX_PIXEL_HEIGHT/8+EXPANSION);
  cairo_move_to(cr,0+EXPANSION+1,5*MAX_PIXEL_HEIGHT/8+EXPANSION);
  cairo_line_to(cr,MAX_PIXEL_WIDTH+EXPANSION-1,5*MAX_PIXEL_HEIGHT/8+EXPANSION);
  cairo_move_to(cr,0+EXPANSION+1,1*MAX_PIXEL_HEIGHT/8+EXPANSION);
  cairo_line_to(cr,MAX_PIXEL_WIDTH+EXPANSION-1,1*MAX_PIXEL_HEIGHT/8+EXPANSION);
  // Lineas en vertical
  cairo_move_to(cr,MAX_PIXEL_HEIGHT+EXPANSION,0+EXPANSION);
  cairo_line_to(cr,MAX_PIXEL_HEIGHT+EXPANSION,MAX_PIXEL_WIDTH+EXPANSION);
  cairo_move_to(cr,MAX_PIXEL_HEIGHT/2+EXPANSION,0+EXPANSION);
  cairo_line_to(cr,MAX_PIXEL_HEIGHT/2+EXPANSION,MAX_PIXEL_WIDTH+EXPANSION);
  cairo_move_to(cr,MAX_PIXEL_HEIGHT/4+EXPANSION,0+EXPANSION);
  cairo_line_to(cr,MAX_PIXEL_HEIGHT/4+EXPANSION,MAX_PIXEL_WIDTH+EXPANSION);
  cairo_move_to(cr,3*MAX_PIXEL_HEIGHT/4+EXPANSION,0+EXPANSION);
  cairo_line_to(cr,3*MAX_PIXEL_HEIGHT/4+EXPANSION,MAX_PIXEL_WIDTH+EXPANSION);
  cairo_move_to(cr,1*MAX_PIXEL_HEIGHT/8+EXPANSION,0+EXPANSION);
  cairo_line_to(cr,1*MAX_PIXEL_HEIGHT/8+EXPANSION,MAX_PIXEL_WIDTH+EXPANSION);
  cairo_move_to(cr,3*MAX_PIXEL_HEIGHT/8+EXPANSION,0+EXPANSION);
  cairo_line_to(cr,3*MAX_PIXEL_HEIGHT/8+EXPANSION,MAX_PIXEL_WIDTH+EXPANSION);
  cairo_move_to(cr,5*MAX_PIXEL_HEIGHT/8+EXPANSION,0+EXPANSION);
  cairo_line_to(cr,5*MAX_PIXEL_HEIGHT/8+EXPANSION,MAX_PIXEL_WIDTH+EXPANSION);
  cairo_move_to(cr,7*MAX_PIXEL_HEIGHT/8+EXPANSION,0+EXPANSION);
  cairo_line_to(cr,7*MAX_PIXEL_HEIGHT/8+EXPANSION,MAX_PIXEL_WIDTH+EXPANSION);
  cairo_stroke(cr);
  cairo_save(cr);
}

int main(int argc, char *argv[])
{
    if (argc!=4)
    {
        cout<<"Argumento faltante..."<<endl;
        exit(1); //Cerramos el programa si es que llega a faltar algun
    }
    double a,b;
    string expresion; //Expresion que contiene la función que deberá ser evaluada
    expresion=argv[1];
    FunctionParser fp;
    fp.AddConstant("pi",3.1415926535897932);
    fp.AddConstant("e",2.718281828459); //Valores más comunes encontrados en matemáticas
    fp.Parse(expresion,"x");
    a=atof(argv[2]);
    b=atof(argv[3]);
    if (a>b) swap(a,b); //Los límites deben de ser de menor a mayor, eso es solo por descuido del usuario o cosas similares
    Info p1(INF,-1,-1);
    Info p2(-INF,-1,-1);
    pair<Info,Info> par=make_pair(p1,p2);
    vector<double> points=utils::linspace(a,b,MAX_POINTS); //The number of points to evaluate the funcion
    vector<double> values=utils::evaluate(points,fp,&par); //Tabién sacamos el máximo y el mínimo de la función y su indice
    //Tenemos ya una evaluación de los puntos en la discretización que necesitábamos por lo que procederemos a graficar
    NZ=utils::ConvertCoordinates(points,values,&par,a,b);
    if(par.second.f<par.first.f) swap(par.second.f,par.first.f);
    Datos.add_values(b,a,par.second.f,par.first.f);
    cairo_t *cr;
    GtkWidget *window;
    GtkWidget *darea;
    gtk_init(&argc, &argv);
    window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    darea = gtk_drawing_area_new();
    gtk_container_add(GTK_CONTAINER(window), darea);
    gtk_window_set_position(GTK_WINDOW(window), GTK_WIN_POS_CENTER);
    gtk_window_set_default_size(GTK_WINDOW(window), MAX_PIXEL_WIDTH+2*EXPANSION, MAX_PIXEL_HEIGHT+2*EXPANSION);
    gtk_window_set_title(GTK_WINDOW(window), "Gráfica de la funcion (Daniel VaAl)");
    g_signal_connect(G_OBJECT(darea),"draw",G_CALLBACK(on_draw_event), cr);
    gtk_widget_show_all(window);
    gtk_main();
    return 0;
}
