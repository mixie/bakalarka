#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <utility>
#include <cstdio>
#include <cassert>
#include <queue>
#include <unordered_set>
#include <string.h>



using namespace std;


#define PI 3.14159265
#define num_colors 255

vector<vector<pair<int,int> > > cir;

/**
Funkcia, ktora vytvori tie "chlieviky" pre rozne vzdialenosti od stredu kruhu
*/
void init_circles(vector<vector<pair<int,int> > >	& cir, int size,int res){
	cir.resize(size*2);
	for(int i=-size;i<size;i++){
		for(int j=-size;j<size;j++){
			int temp=0.5+sqrt(j*j+i*i);
			cir[temp].push_back( make_pair(i,j));
		}
	}
	for(int i=0;i<size;i++){
		random_shuffle(cir[i].begin(),cir[i].end());
		if(cir[i].size()>res){
			cir[i].resize(res);
		}

	}
}


struct picture{
	int x,y;
	int rad=-1;
	vector<vector<double> > m;
	picture();
	picture(char *);
	picture(int s_x, int s_y);
	picture(int s_x,int s_y,double def);
	void save(char *); 
	bool inrange(int, int);
	void setRad(int);
};

picture::picture(char * file){
	FILE *fr=fopen(file,"r");
	char str[500];
	fgets (  str, 500, fr);
	fgets ( str, 500, fr);
	int colors;
	fscanf(fr,"%d %d",&x,&y);
	m.resize(y+47);
	for(int i=0;i<y;i++){
		for(int j=0;j<x;j++){
			double num;
			fscanf(fr,"%lf",&num);
			m[i].push_back(num/num_colors);
		}
	}
	fclose(fr);
}


picture::picture(int s_x,int s_y,double d){
	x=s_x;
	y=s_y;	
	m.resize(y+47);
	for(int i=0;i<y;i++){
		for(int j=0;j<x;j++){
			m[i].push_back(d);
		}
	}
}


picture::picture(int s_x, int s_y):picture(s_x,s_y,0.0) {}

picture::picture(){

}

void picture::setRad(int r){
	rad=r;
}


void picture::save(char * file){
	FILE *fw=fopen(file,"w");
	fprintf(fw, "P2\n");
	fprintf(fw, "%d %d 255\n",x,y);	
	//printf("%f\n", m.size());
	for(int i=0;i<y;i++){
		for(int j=0;j<x;j++){
			int a=(int)255.0*m[i][j];
			if(a<0){
				fprintf(fw, "0\n");
			}else{
				fprintf(fw, "%d\n", a);
			}
		}
	}
	fclose(fw);
}

bool picture::inrange(int p_x, int p_y){
	if((p_x>=0)&&(p_y>=0)&&(p_x<x)&&(p_y<y)){
		return true;
	}
	return false;
}


struct cilia{
	double m=-1.;
	picture * pic;
	int rad;
	int cen_x,cen_y;
	cilia(int, int, picture * );
	int ret_xy(int , int );
	double simple_repres(vector<double> & );
	double mean();
};

cilia::cilia(int x, int y,picture * p){
	cen_x=x;
	cen_y=y;
	pic=p;
	rad=p->rad;
}


int cilia::ret_xy(int x, int y){
	return pic->m[cen_y-rad+y][cen_x-rad+x];
}

double cilia::simple_repres(vector<double> & simple_repres){
	simple_repres.resize(rad*2+47,0.);
	for(int i=0;i<rad;i++){
		for(auto p:cir[i]){
			int y, x; 
			tie(x,y) = p;
			simple_repres[i] += pic->m[cen_y+y][cen_x+x];
		}
		simple_repres[i]=simple_repres[i]/cir[i].size();
	}
	return accumulate(simple_repres.begin(),simple_repres.end(),0.)/(rad);
}


double cilia::mean(){
	if(m==-1){
		double sum=0.;
		for(int i=0;i<rad*2;i++){
			for(int j=0;j<rad*2;j++){
				sum+=ret_xy(j,i);
			}
		}
		m=sum/((rad*2)*(rad*2));
	}return m;
}

double simple_pearson(cilia & a, cilia & b){
	vector<double> simple_repres_a,simple_repres_b;
	double mean_a=a.simple_repres(simple_repres_a);
	double mean_b=b.simple_repres(simple_repres_b);
	assert(a.rad==b.rad);
	double sum=0; double d_a=0; double d_b=0;
	for(int i=0;i<a.rad;i++){
		sum+=(simple_repres_a[i]-mean_a)*(simple_repres_b[i]-mean_b);
		d_a+=(simple_repres_a[i]-mean_a)*(simple_repres_a[i]-mean_a);
		d_b+=(simple_repres_b[i]-mean_b)*(simple_repres_b[i]-mean_b);
	}
	return sum/(sqrt(d_a)*sqrt(d_b));
}

double pearson(cilia & a, cilia & b){
	double mean_a=a.mean();
	double mean_b=b.mean();
	assert(a.rad==b.rad);
	double sum=0; double d_a=0; double d_b=0;
	for(int i=0;i<2*a.rad;i++){
		for(int j=0;j<2*a.rad;j++){
			sum+=(a.ret_xy(j,i)-mean_a)*(b.ret_xy(j,i)-mean_b);
			d_a+=(a.ret_xy(j,i)-mean_a)*(a.ret_xy(j,i)-mean_a);
			d_b+=(b.ret_xy(j,i)-mean_b)*(b.ret_xy(j,i)-mean_b);
		}
	}

	return sum/(sqrt(d_a)*sqrt(d_b));
}

//TODO: treba tam tych 10?
bool inrange(int x, int y, cilia & c, picture * p){
	if((x>p->x-c.rad)||(x<c.rad+10)||(y>p->y-c.rad)||(y<c.rad+10)){
		return false;
	}
	return true;
}


/**
Funkcia, ktora pre jednu sample riasinku skusi najst ine riasinky na obrazku
sam - sample riasinka
in,out - vstupny, vystupny obrazok
int simple - (0/1) 0 - rata korelaciu pre celu riasinku, 1 - rata korelaciu sum bodov podla vzdialenosti od stredu 
*/
void pear(cilia sam,picture * in,picture * out,int simple){
	for(int i=0;i<in->y;i++){
		printf("%d\n", i);
		for(int j=0;j<in->x;j++){
			if(inrange(j,i,sam,in)){
				cilia act(j,i,in);
				if(simple){
					out->m[i][j]=simple_pearson(sam,act);
				}else{
					out->m[i][j]=pearson(sam,act);	
				}
			}else
			out->m[i][j]=0;
		}
	}
}

void pear_selective(cilia sam,picture * in1,picture * in2, picture * out,int simple, double threshold){
	vector <pair<double,pair<int,int > > >to_sort; 
	for(int i=0;i<in2->y;i++){
		for(int j=0;j<in2->x;j++){
			if(in2->m[i][j]<threshold){
				to_sort.push_back(make_pair(in2->m[i][j],make_pair(i,j)));
			}
		}
	}
	sort(to_sort.begin(),to_sort.end());
	printf("COUNTER\n");
	int counter=0;
	for(int i=0;i<to_sort.size();i++){
		counter++;
		if(counter%10000==0){
			printf("%d\n", counter);
		}
		if(to_sort[i].first<threshold){
			if(inrange(to_sort[i].second.second,to_sort[i].second.first,sam,in1)){
				cilia act(to_sort[i].second.second,to_sort[i].second.first,in1);
				//printf("%d %d\n",to_sort[i].second.second,to_sort[i].second.first);
				out->m[to_sort[i].second.first][to_sort[i].second.second]=simple_pearson(sam,act);
			}
		}else{
			return;
		}
	}
}


bool isNewCiliaCentre(int x1,int y1,int x2,int y2, int rad){
//	printf("%d %d %d %d %f\n", x1,x2,y1,y2,sqrt(abs(x1-x2)*abs(x1-x2)+abs(y1-y2)*abs(y1-y2)));
	return sqrt(abs(x1-x2)*abs(x1-x2)+abs(y1-y2)*abs(y1-y2))>2*rad;
}

struct pair_hash {
	inline std::size_t operator()(const std::pair<int,int> & v) const {
		return v.first*31+v.second;
	}
};
/**
Ak je najdeny priblizny stred riasinky, tato funkcia spusti BFSko 
na jeho okolie, ktore je dostatocne biele (t.j. dostatocne stred) a najde priemer
in - vstupnyobrazok
x,y - suradnice zatial najdeneho stredu
threshold - hranica, po ktoru povazuje body za mozne stredy
*/
pair<int,int> findExactCentre(picture * in, int x, int y, double threshold){
	queue <pair<int,int> > q;
	unordered_set <pair<int,int>,pair_hash > visited;
	q.push(make_pair(x,y));
	visited.insert(make_pair(x,y));
	int dx[]={1,-1,0,0};
	int dy[]={0,0,1,-1};
	int sumX=0,sumY=0;
	int count=0;
	while(!q.empty()){
		tie(x,y)=q.front();
		sumX+=x;
		sumY+=y;
		count++;
		q.pop();
		for(int i=0;i<4;i++){
			if(in->inrange(x+dx[i],y+dy[i])){
				if((in->m[y+dy[i]][x+dx[i]]>threshold)&&(visited.count(make_pair(x+dx[i],y+dy[i]))==0)){
					q.push(make_pair(x+dx[i],y+dy[i]));
					visited.insert(make_pair(x+dx[i],y+dy[i]));
				}
			}
		}
	}
	return make_pair(sumY/count,sumX/count);
}

/**
Nakresli trosku vacsi bod do obrazku
*/
void draw_point(int rad,int x, int y, picture * p,double color){
	for(int i=y-rad;i<y+rad;i++){
		for(int j=x-rad;j<x+rad;j++){
			if(p->inrange(j,i)){
				p->m[i][j]=color;
			}
		}
	}
}



/**
 Hlada centra riasiniek z predspracovaneho orazku, zoradi si body podla "belosti" a potom hlada take, co su dostatocne
 daleko od seba, aby mohli byt stredy a tie este potom vylepsuje cez findExactCircles
 in - vstupny obrazok, ktory bol ako vystup funkcie pear
out -vystupny obrazok, kdesu iba nakreslene bodky na mieste riasiniek
rad- polomer riasinky 
threshold - po aku hranicu este bude povazovat body za mozny stred
centres - vysledny vektor stredov riasiniek
*/
void findCentres(picture * in, picture * out, int rad,double threshold, vector<pair<int,int>> & centres){
	threshold=threshold/255.0;
	vector <pair<double,pair<int,int > > >to_sort; 
	for(int i=0;i<in->y;i++){
		for(int j=0;j<in->x;j++){
			if(in->m[i][j]!=0.0){
				to_sort.push_back(make_pair(in->m[i][j],make_pair(i,j)));
			}
		}
	}
	sort(to_sort.begin(),to_sort.end());
	for(int i=to_sort.size()-1;i>=0;i--){
		int an=true;
		for(int j=0;j<centres.size();j++){
			if(!isNewCiliaCentre(to_sort[i].second.first,to_sort[i].second.second,
				centres[j].first,centres[j].second,rad)){
				an=false;
				break;
			}
		}
		if(an){
			if(to_sort[i].first>threshold){
				//printf("%d %d %f\n", to_sort[i].second.second,to_sort[i].second.first,to_sort[i].first);
				pair<int, int> p=findExactCentre(in,to_sort[i].second.second,to_sort[i].second.first,threshold);
				centres.push_back(p);
			}
		}
	}
	for(int i=0;i<centres.size();i++){
		draw_point(3,centres[i].second,centres[i].first,out,0.0);
	}
}
/**
Najde orientaciu riasinky, na vstupe potrebuje povodny obrazok a presny stred kazdej riasinky,
zatial to len vykresluje
in - vstupny obrazok, mal by to byt povodny, nespracovany obrazok riasinky
out -vystupny obrazok
rad- polomer riasinky 
circles - vysledny vektor stredov riasiniek
num_lines - koľko veľa rôzne otočených priamok bude skúšať
**/
void findOrientation(picture * in, picture * out, int rad, vector<pair<int,int> > & centres,int num_lines){
	rad=rad/4;
	out->m.resize(out->y);
	for(int i=0;i<out->y;i++){
		for(int j=0;j<out->x;j++){
			out->m[i].push_back(1.0);
		}
	}
	for(auto &v:centres){
		int x,y;
		tie(y,x)=v;
		double min=99999999999;
		double m_x=0,m_y=0;
		for(int k=0;k<num_lines;k++){
			double sum=0;
			pair<double,double> ab,ac;
			ab.second=cos(((PI)/num_lines)*k);
			ab.first=sin(((PI)/num_lines)*k);
			for(int i=y-rad;i<y+rad;i++){
				for(int j=x-rad;j<x+rad;j++){
					if(sqrt((x-j)*(x-j)+(y-i)*(y-i))<rad){
						ac.first=i-y; ac.second=j-x;
						if(in->inrange(j,i)){
							sum+=(in->m[i][j]*255)/(abs(-ab.second*ac.first+ab.first*ac.second)+1);
						}
					}
				}
			}
			if(sum<min){
				min=sum;
				m_x=ab.first;
				m_y=ab.second;
			}
		}
		for(int i=-25;i<25;i+=2){
			draw_point(3,x+m_y*i,y+m_x*i,out,0.0);
		}
	}
}

void preprocessPrefix(picture * in, picture * out){
	vector <vector<double > > pr;
	pr.resize(in->y+5);
	for(int i =0;i<in->y+5;i++){
		pr[i].resize(in->x+5);
	}
	printf("%d %d\n", in->y,in->x);
	for(int i=0;i<in->y;i++){
		pr[i][0]=0;
	}
	for(int j=0;j<in->x;j++){
		pr[0][j]=0;
	}
	for(int i=1;i<=in->y;i++){
		for(int j=1;j<=in->x;j++){
			pr[i][j]=pr[i][j-1]+pr[i-1][j]-pr[i-1][j-1]+in->m[i-1][j-1];
		}
	}
	int rad=in->rad;
	double max=0;
	for(int i=rad;i<in->y-rad;i++){
		for(int j=rad;j<in->x-rad;j++){
			out->m[i][j]=pr[i+1+rad][j+1+rad]-pr[i+1+rad][j+1-rad]-pr[i+1-rad][j+1+rad]+pr[i+1-rad][j+1-rad];
			if(out->m[i][j]>max){
				max=out->m[i][j];
			}
		}
	}
	for(int i=rad;i<in->y-rad;i++){
		for(int j=rad;j<in->x-rad;j++){
			out->m[i][j]=out->m[i][j]/max;
		}
	}
}

int threshold1(picture * in, picture * out, double threshold){
	int c=0;
	for(int i=0;i<in->y;i++){
		for(int j=0;j<in->x;j++){
			if(in->m[i][j]<threshold){
				c++;
			}
			out->m[i][j]=(int)(in->m[i][j]>threshold);
		}
	}
	return c;
}

int main(int argc, char * argv[]){
/**	int cilia_rad=80;
	int random_points=20;
	int threshold=200;
	int num_lines=100;
	init_circles(cir,cilia_rad,random_points);
	picture p1("rias.pgm");
	p1.setRad(cilia_rad);
	picture p2(p1.x,p1.y,0.0);
	cilia c(346,156,p1); //x=135, y=100
	pear(c,p1,p2,1);
	p2.save("rias-out-pom1.pgm");
	picture p3(p1.x,p1.y,1.0);
	vector <pair<int,int>> centres;
	findCentres(p2,p3,cilia_rad,threshold,centres);
	p3.save("rias-out-pom2.pgm");
	picture p4(p1.x,p1.y,1.0);
	findOrientation(p1,p4,cilia_rad,centres,num_lines);
	p4.save("rias-out.pgm");
	vector <pair<int,int>> centres2;
	picture p6("Tv4-small-out-crop.pgm");
	picture p5("Tv4-small-crop.pgm");
	p5.setRad(80);
	p6.setRad(80);
	picture p7(p5.x,p5.y,1.0);
	picture p8(p5.x,p5.y,1.0);
	findCentres(p6,p7,cilia_rad,threshold,centres2);
	findOrientation(p5,p8,cilia_rad,centres2,num_lines);
	p7.save("Tv4-small-crop-centres.pgm");
	p8.save("Tv4-small-crop-orient.pgm");**/
/*	int cilia_rad=80;
	int random_points=10;
	int threshold=200;
	int num_lines=100;
	init_circles(cir,cilia_rad,random_points);
	picture p1("Tv4-small-crop.pgm");
	p1.setRad(cilia_rad);
	picture p2(p1.x,p1.y,1);
	printf("%d %d\n", p1.x,p1.y);
	preprocessPrefix(p1,p2);
	p2.save("rias-out.pgm");
	picture p4(p1.x,p1.y,1);
	printf("%d",threshold1(p2,p4,0.62));
	p4.save("threshold.pgm");
 	picture p3(p1.x,p1.y,0);
	cilia c(267,175,p1);
	pear_selective(c,p1,p2, p3,1, 0.62);
	p3.save("rias-out3.pgm");
	picture p5(p1.x,p1.y,1.0);
	vector <pair<int,int>> centres;
	findCentres(p3,p5,cilia_rad,threshold,centres);
	p5.save("rias-out-pom23.pgm");
	picture p6(p1.x,p1.y,1.0);
	findOrientation(p1,p6,cilia_rad,centres,num_lines);
	p6.save("rias-out23.pgm");*/
	int cilia_rad=80;
	int random_points=10;
	int threshold=200;
	int num_lines=100;
	int x=1377,y=1762;
	double threshold2=0.52;
	init_circles(cir,cilia_rad,random_points);
	char s[100]="fotky/Tv4_1"; //fotky/Tv10 x=1814, y=1230
	char s2[100];
	strcpy(s2,s);
	picture * p1= new picture(strcat(s,".pgm"));
	strcpy(s,s2);
	p1->setRad(cilia_rad);
	picture * p2= new picture(p1->x,p1->y,1);
	preprocessPrefix(p1,p2);
	p2->save(strcat(s,"02.pgm"));
strcpy(s,s2);	
picture * p3= new picture(p1->x,p1->y,1);
	printf("points:%d\n",threshold1(p2,p3,threshold2));
	p3->save(strcat(s,"03.pgm"));
strcpy(s,s2);	
cilia c(x,y,p1);
	picture * p4= new picture(p1->x,p1->y,0);
	pear_selective(c,p1,p2,p4,1,threshold2);
	p4->save(strcat(s,"04.pgm"));
strcpy(s,s2);	
picture * p5= new picture(p1->x,p1->y,1);
	vector <pair<int,int>> centres;
	findCentres(p4,p5,cilia_rad,threshold,centres);
	p5->save(strcat(s,"05.pgm"));
strcpy(s,s2);	
picture * p6= new picture(p1->x,p1->y,1.0);
	findOrientation(p1,p6,cilia_rad,centres,num_lines);
	p6->save(strcat(s,"06.pgm"));
	delete(p1);	delete(p2);	delete(p3);
delete(p4);	delete(p5);	delete(p6);

}
