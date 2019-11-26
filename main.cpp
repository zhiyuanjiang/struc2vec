#include <iostream>
#include <thread>
#include <queue>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include "FastDTW.h"
#include "Dist.h"

using namespace std;
using namespace fastdtw;

//节点个数
const int N = 10000;
//线程个数
const int M = 4;
//"邻居"个数
const int NUM_NEIGHBOR = 2.0*log(N)/log(2);

const int INF = 0x3ffffff;

//边
struct Edge{
    int u,v;
    double val;
    Edge(int _u, int _v, double _val):u(_u),v(_v),val(_val){}
};

//节点的度
struct Degree{
    int d, id;
    bool operator < (const Degree&a) const{
        return d < a.d;
    }
}degree[N];

//a[0]:节点的度d，a[1]:节点的度为d的个数
struct Node{
    double a[2];
    Node(int v1, int v2){
        a[0] = v1; a[1] = v2;
    }
};


//原始图
//vector<vector<int>>graph = {{1,2},{0,3,4},{0},{1},{1}};
vector<int>graph[N];
//节点某一阶的向量
//vector<vector<vector<int>>>simi_graph;
vector<vector<vector<Node>>>simi_graph;
//节点u的邻居,u只需要和他的“邻居”计算,不需要和所有的节点计算
vector<int>neighbors[N];
//存储生成图的边
vector<Edge>edge[N/M+1];
//多线程
int pool[N/M+1];

int vis[N], L[N];

bool read_file(string& filename)
{
    ifstream file(filename.c_str(), ifstream::in);

    if(!file){
        cout << "failed to open file "+filename+"!" << endl;
        return false;
    }

    string line;
    int u,v;
    while(!file.eof()){
        getline(file, line);
        istringstream record(line);
        record >> u >> v;
        //cout << u << " " << v << endl;
        graph[u].push_back(v);
        graph[v].push_back(u);
    }

    file.close();
    return true;
}

//O(nlogn)
void get_neighbors()
{
    for(int i = 0; i < N; ++i){
        degree[i].id = i;
        degree[i].d = graph[i].size();
    }
    sort(degree, degree+N);

//    for(int i = 0; i < N; ++i){
//        cout << degree[i].d << " ";
//        if(i%10 == 0) cout << endl;
//    }

    for(int i = 0; i < N; ++i){
        int d = graph[i].size();
        //int p = lower_bound()

        int l = 0, r = N;
        int mid, ans;
        while(l <= r){
            mid = (l+r)>>1;
            if(degree[mid].d <= d){
                ans = mid;
                l = mid+1;
            }else{
                r = mid-1;
            }
        }

        int cnt = 0, a, b;
        l = ans-1, r = ans+1;
        if(degree[ans].id != i){
            neighbors[i].push_back(degree[ans].id);
            cnt++;
        }
        while(cnt < NUM_NEIGHBOR){
            a = (a<0)  ? INF : degree[l].d;
            b = (b>=N) ? INF : degree[r].d;
            if(abs(a-d) < abs(b-d)){
                neighbors[i].push_back(degree[l].id);
                l--;
            }else{
                neighbors[i].push_back(degree[r].id);
                r++;
            }
            cnt++;
        }
    }
}

bool write_file(string& filename, vector<Edge>& edge)
{
    ofstream file(filename.c_str(), ofstream::app);

    if(!file){
        cout << "failed to open file "+filename+"!" << endl;
        return false;
    }

    ostringstream data;
    for(auto it:edge){
        string v = to_string(it.val);
        v = v.substr(0, v.length()-4);
        string s = to_string(it.u)+" "+to_string(it.v)+" "+v;
        //cout << s << endl;
        data << s << endl;
        //data << it.u << " " << it.v << " " << to_string(it.val) << endl;
        //cout << to_string(it.val) << endl;
    }

    file.write(data.str().c_str(), data.str().length());

    file.close();
    return true;
}

void change(vector<int>&a, vector<Node>&b)
{
    int cnt = 0,flag = a[0];
    for(int i = 0; i < a.size(); ++i){
        if(flag == a[i]) cnt++;
        else{
            b.push_back(Node(flag, cnt));
            cnt = 1; flag = a[i];
        }
    }
    b.push_back(Node(flag, cnt));
}

//O(n)
vector<vector<Node>> bfs(int root, int pos)
{
    vector<vector<Node>> a;
    vector<int> b;

    vector<Node> c;

    memset(vis, 0, sizeof(vis));
    memset(L, 0, sizeof(L));
    queue<int> q;
    q.push(root);
    vis[root] = 1;
    int tick = 0;
    while(!q.empty()){
        int u = q.front();
        q.pop();
        if(L[u] == tick){
            b.push_back(graph[u].size());
        }else{
            sort(b.begin(), b.end());
            change(b, c);
            a.push_back(c);
            b.clear();

            c.clear();

            tick = L[u];
            b.push_back(graph[u].size());
        }
        if(L[u] == pos-1) continue;
        for(int i = 0; i < graph[u].size(); ++i){
            int v = graph[u][i];
            if(vis[v]) continue;
            q.push(v);
            vis[v] = 1;
            L[v] = L[u]+1;
        }
    }
    sort(b.begin(), b.end());
    change(b, c);
    a.push_back(c);
    return a;
}

double dtw(vector<double>&a, vector<double>&b)
{
    int r = a.size(), c= b.size();
    double dp[r+1][c+1];
    memset(dp, 0, sizeof(dp));
    for(int i = 1; i <= r; ++i){
        for(int j = 1; j <= c; ++j){
            dp[i][j] = max(a[i-1],b[j-1])*1.0/min(a[i-1],b[j-1])-1;
            dp[i][j] = min(min(dp[i-1][j-1]+2*dp[i][j],dp[i-1][j]+dp[i][j]),min(dp[i-1][j]+dp[i][j],dp[i][j-1]+dp[i][j]));
        }
    }
    return dp[r][c];
}

double dist(Node& a, Node& b)
{
    double exp = 0.5;
    double mi = min(a.a[0], b.a[0])+0.5;
    double ma = max(a.a[0], b.a[0])+0.5;
    return ma/mi-1;
}

double myfastdtw(vector<Node>&a, vector<Node>&b)
{
    //double *p1 = &a[0], *p2 = &b[0];
    TimeSeries<double,2> tsI;
    for(int i = 0; i < a.size(); ++i) {
        tsI.addLast(i, TimeSeriesPoint<double,2>(a[i].a));
    }

    TimeSeries<double,2> tsJ;
    for(int i = 0; i < b.size(); ++i)
    {
        tsJ.addLast(i, TimeSeriesPoint<double,2>(b[i].a));
    }

    double ans = 0;
    if(a.size() == 1 && b.size() == 1) ans = dist(a[0], b[0]);
    else if(a.size() < 100 || b.size() < 100) ans = STRI::getWarpDistBetween(tsI, tsJ, Dist());
    else ans = FAST::getWarpDistBetween(tsI, tsJ, JInt(1), Dist());

    return ans;
    //return STRI::getWarpDistBetween(tsI, tsJ, Dist());
    //return FAST::getWarpDistBetween(tsI, tsJ, JInt(1), Dist());
}

//O(nlogn*l)
void buildGraph(int s, int e, int pos, int index)
{
//    vector<Node>a, b;
    for(int i = s; i < e; ++i){
        for(auto j : neighbors[i]){
            if(i == j || simi_graph[i].size() < pos || simi_graph[j].size() < pos) continue;

//            a.clear(); b.clear();
//            for(auto it : simi_graph[i][pos]) a.push_back(it);
//            for(auto it : simi_graph[j][pos]) b.push_back(it);

            double val = myfastdtw(simi_graph[i][pos], simi_graph[j][pos]);

//            cout << i << " " << j << " " << a.size() << " " << b.size() << endl;
//            for(auto it : a) cout << it.a[0] << " " << it.a[1] << endl;
//            cout << endl;
//            for(auto it : b) cout << it.a[0] << " " << it.a[1] << endl;;
//            cout << endl;
//            cout << val << endl;
//            return;

            edge[index].push_back(Edge(i, j, val));
        }
    }
    pool[index] = 1;
}

//代码应该还是存在问题，还要要看懂fastdtw，fastdtw可能没有使用正确。

int main()
{
    string graph_file = "f://test2.txt";
    string output = "f://";

    cout << "test" << endl;
    time_t s,e;

    s = time(0);

    //读入图
    read_file(graph_file);
    //获取每个节点的邻居
    get_neighbors();

    for(int i = 0; i < N; ++i){
        simi_graph.push_back(bfs(i, 3));
    }

    e = time(0);

    cout << "deal data:" << (e-s) << endl;

    //buildGraph(0, N-1, 2, 0);
    //return 0;

    for(int x = 0; x < 3; ++x){

        s = time(0);

        memset(pool, 0, sizeof(pool));
        int cnt = 0;
        for(int i = 0; i < N; i = i+N/M, ++cnt){
            if(i+N/M >= N-1){
                thread t(buildGraph, i, N-1, x, cnt);
                t.detach();
                //pool.push_back(thread(buildGraph, i, N-1, x, j));
            }else{
                thread t(buildGraph, i, i+N/M, x, cnt);
                t.detach();
                //pool.push_back(thread(buildGraph, i, i+N/M, x, j));
            }
        }

        while(1){
            int flag = 1;
            for(int i = 0; i < cnt; ++i)
                flag &= pool[i];
            if(flag) break;
            //else std::this_thread::sleep_for(chrono::seconds(1));
        }

        e = time(0);
        cout << "get the edge:" << (e-s) << endl;

        s = time(0);

        for(int i = 0; i < cnt; ++i){
            string name = output+"edgelist"+to_string(x)+".txt";
            write_file(name, edge[i]);
            edge[i].clear();
        }

        e = time(0);
        cout << "write the edge:" << (e-s) << endl;

    }

    /*
    100 5 24 70 0
    1000 5 80 600 5
    10000 5 100 2500 40
    100000 5 100 2500 530
    */

    return 0;
}
