#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <utility>
#include <cstdio>
#include <cassert>
#include <queue>
#include <unordered_set>


using namespace std;


vector<vector<pair<int,int> > > cir;

#define PI 3.14159265

void init_circles(vector<vector<pair<int,int> > >	& cir, int size){
	cir.resize(2*size);
	for(int i=-size;i<size;i++){
		for(int j=-size;j<size;j++){
			int temp=0.5+sqrt(j*j+i*i);
			cir[temp].push_back( make_pair(i,j));
		}
	}
}

struct picture{
	int x,y;
	vector<vector<double> > m;
	picture();
	picture(char *);
	picture(int s_x, int s_y);
	void save(char *);
	bool inrange(int, int);
	inrange(int x, int y, cilia & c);
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
			int num;
			fscanf(fr,"%d",&num);
			m[i].push_back(num);
		}
	}
	fclose(fr);
}

picture::picture(int s_x, int s_y){
	x=s_x;
	y=s_y;	
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

picture::picture(){

}

struct cilia{
	double m=-1.;
	picture pic;
	int rad;
	int cen_x,cen_y;
	cilia(int, int, int ,picture & );
	int ret_xy(int , int );
	double simple_repres(vector<double> & );
	double mean();
};

cilia::cilia(int x, int y, int s,picture & p){
	cen_x=x;
	cen_y=y;
	rad=s;
	pic=p;
}

int cilia::ret_xy(int x, int y){
	assert(cen_y-rad+y>=0);
	assert(cen_x-rad+x>=0);
	assert(cen_y-rad+y<pic.y);
	assert(cen_x-rad+x<pic.x);
	return pic.m[cen_y-rad+y][cen_x-rad+x];
}

double cilia::simple_repres(vector<double> & simple_repres){
	simple_repres.resize(rad*2+47,0.);
	for(int i=0;i<rad;i++){
		for(auto p:cir[i]){
			int y, x; 
			tie(x,y) = p;
			simple_repres[i] += pic.m[cen_y+y][cen_x+x];
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

picture::inrange(int x, int y, cilia & c, picture & p){
	if((x>p.x-c.rad)||(x<c.rad+10)||(y>p.y-c.rad)||(y<c.rad+10)){
		return false;
	}
	return true;
}

void pear(cilia sam,picture & in,picture & out,int simple){
	out.m.resize(in.m.size());
	for(int i=0;i<in.y;i++){
		printf("%d\n", i);
		for(int j=0;j<in.x;j++){
			if(inrange(j,i,sam,in)){
				cilia act(j,i,sam.rad,in);
				if(simple){
					out.m[i].push_back(simple_pearson(sam,act));
				}else{
					out.m[i].push_back(pearson(sam,act));	
				}
			}else{
				out.m[i].push_back(0);
			}
		}
	}
}

bool isNewCiliaCentre(int x1,int y1,int x2,int y2, int rad){
	return sqrt(abs(x1-x2)*abs(x1-x2)+abs(y1-y2)*abs(y1-y2))>2*rad;
}

struct pair_hash {
	inline std::size_t operator()(const std::pair<int,int> & v) const {
		return v.first*31+v.second;
	}
};

pair<int,int> findExactCentre(picture & in, int x, int y,int rad, int threshold){
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
			if(in.inrange(x+dx[i],y+dy[i])){
				if((in.m[y+dy[i]][x+dx[i]]>threshold)&&(visited.count(make_pair(x+dx[i],y+dy[i]))==0)){
					q.push(make_pair(x+dx[i],y+dy[i]));
					visited.insert(make_pair(x+dx[i],y+dy[i]));
				}
			}
		}
	}
	return make_pair(sumY/count,sumX/count);
}

void draw_point(int rad,int x, int y, picture & p,double color){
	for(int i=y-rad;i<y+rad;i++){
		for(int j=x-rad;j<x+rad;j++){
			if((0<=i)&&(0<=j)&&(i<p.y)&&(j<p.x)){
				p.m[i][j]=color;
			}
		}
	}
}


void findCentres(picture & in, picture & out, int rad,int threshold, vector<pair<int,int>> & centres){
	vector <pair<double,pair<int,int > > >to_sort; 
	for(int i=0;i<in.y;i++){
		for(int j=0;j<in.x;j++){
			if(in.m[i][j]!=0.0){
				to_sort.push_back(make_pair(in.m[i][j],make_pair(i,j)));
			}
		}
	}
	vector<pair<int,int> > centres2;
	sort(to_sort.begin(),to_sort.end());
	for(int i=to_sort.size()-1;i>=0;i--){
		int an=true;
		for(int j=0;j<centres2.size();j++){
			if(!isNewCiliaCentre(to_sort[i].second.first,to_sort[i].second.second,
				centres2[j].first,centres2[j].second,80)){
				an=false;
			break;
		}
	}
	if(an){
		if(to_sort[i].first>threshold){
			pair<int, int> p=findExactCentre(in,to_sort[i].second.second,to_sort[i].second.first,80,threshold);
			centres2.push_back(p);
		}
	}
}
out.m.resize(out.y);
for(int i=0;i<out.y;i++){
	for(int j=0;j<out.x;j++){
		out.m[i].push_back(1.0);
	}
}
for(int i=0;i<centres2.size();i++){
	draw_point(3,centres2[i].second,centres2[i].first,out,0.0);
}
centres=centres2;
}

void lines(picture & in, picture & out, int rad, vector<pair<int,int> > & circles,int num_lines){
	out.m.resize(out.y);
	for(int i=0;i<out.y;i++){
		for(int j=0;j<out.x;j++){
			out.m[i].push_back(1.0);
		}
	}
	for(auto &v:circles){
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
						if((i>=0)&&(j>=0)&&(in.y>i)&&(in.x>j)&&(i!=y)&&(j!=x)){
							double temp=abs(-ab.second*ac.first+ab.first*ac.second);
							//if(temp<0.5){
						//	printf("%f %f %f %f %f %f %f\n",temp,ac.first,ac.second,ab.first,ab.second,in.m[i][j],in.m[i][j]);
							//	sum+=in.m[i][j];
							//}else{
								sum+=in.m[i][j]/(temp+1);
						//	printf("%f %f %f %f %f %f %f\n",temp,ac.first,ac.second,ab.first,ab.second,in.m[i][j],in.m[i][j]/temp);

							//}
						}else{
						//printf("skaasdkdaskASSSSSSSSSSSSSSSSSS\n");
						}
					}
				}
			}
			if(sum<min){
				min=sum;
				m_x=ab.first;
				m_y=ab.second;
			}
		//printf("sum:%f max:%f\n", sum,max);
			printf("baf:  %f %f %f\n", ab.first,ab.second,sum);

		}
		printf("baf:***  %f %f %f %d %d\n", m_x,m_y,min,x,y);
		for(int i=-25;i<25;i+=2){
			draw_point(3,x+m_y*i,y+m_x*i,out,0.0);
		}
	}
}

int main(int argc, char * argv[]){
/**	init_circles(cir,80);
	picture p1("fotky/test-obratene2.pgm");
	/**cilia c(336,113,80,p1);
	picture p2(p1.x,p1.y);
	pear(c,p1,p2,true);
	printf("step1\n");
	picture p3(p1.x,p1.y);
	vector<pair<int, int>> centres;
	findCentres(p2,p3,80,180,centres);
	printf("step2\n");**/
//	vector<pair<int,int>> centres;
	/**centres.push_back(make_pair(392,392));
	centres.push_back(make_pair(728,728));
	centres.push_back(make_pair(932,200));
	centres.push_back(make_pair(92,222));**/
//centres.push_back(make_pair(148,341));

//	picture p4(p1.x,p1.y);
//	lines(p1,p4,20,centres,100);
//	printf("step3\n");
//	p4.save("fotky/lines-obr2.pgm");
	/**picture p1("fotky/Tv9-small-out2.pgm");
	picture p2(p1.x,p1.y);
	vector<pair<int, int>> centres;
	findCentres(p1,p2,80,180,centres);
	p2.save("fotky/Tv9-small-out2-centres.pgm");
	picture p3("fotky/Tv9-small2.pgm");
	picture p4(p3.x,p3.y);
	lines(p3,p4,20,centres,50);
	p4.save("fotky/lines2-min-Tv9.pgm");*/
	/**if(argc!=7){

		printf("%d",argc);
		return 1;
	}
	char ** n;
	picture p1(argv[1]);
	picture p2(p1.x,p1.y);
	int r_x,r_y,rad,an;
	sscanf(argv[2],"%d",&rad);
	sscanf(argv[3],"%d",&r_x);
	sscanf(argv[4],"%d",&r_y);
	sscanf(argv[5],"%d",&an);
	init_circles(cir,rad);
	cilia c(r_x, r_y,rad ,p1 );
//	cilia d(318, 287, 80 ,p1 );
	pear(c,p1,p2,an);
	p2.save(argv[6]);
//	printf("%f \n", simple_pearson(c,d));**/

	picture p1("fotky/Tv4-small-crop.pgm");
	picture p3("fotky/Tv4-small-out-crop.pgm");
	picture p2(p3.x,p3.y);
	picture p4(p1.x,p1.y);
	vector<pair<int,int>> centres;
	//centres.push_back(make_pair(676,400));
	findCentres(p3,p2,80,180,centres);
	lines(p1,p4,20,centres,80);
	p4.save("fotky/Tv4-small-orient-crop.pgm");

}

//tv4small - 543 855
//tv14small - 636 394
//Tv4small3 - 244 520

