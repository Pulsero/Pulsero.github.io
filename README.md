```
// P1036 [NOIP 2002 普及组] 选数
// 已知 n 个整数 x  1 ​  ,x  2 ​  ,⋯,x  n ​  ，以及 1 个整数 k（k<n）。从 n 个整数中任选 k 个整数相加，可分别得到一系列的和。例如当 n=4，k=3，4 个整数分别为 3,7,12,19 时，可得全部的组合与它们的和为
#include <bits/stdc++.h>

inline bool isprime(int x) { // 判断一个数是否是素数
	if (x == 1) return false; // 注意这步特判是必需的
	for (int i = 2; i * i <= x; ++i)
		if (x % i == 0)
			return false;
	return true;
}

const int N = 25;
int a[N], ans, n, k;

void dfs(int now, int sum, int sid) {
	// 现在已经选了 now 个数，当前总和为 sum
	// sid 是这次选数的起始下标，即我们从 a[sid] 开始选数枚举
	if (now == k) {
		if (isprime(sum))
			++ans;
		return ;
	}
	for (int i = sid; i <= n - k + now + 1; ++i)
		dfs(now + 1, sum + a[i], i + 1);
	return ;
}

int main() {
	scanf("%d%d", &n, &k);
	for (int i = 1; i <= n; ++i)
		scanf("%d", &a[i]);
	dfs(0, 0, 1);
	printf("%d\n", ans);
	return 0;
}



// P1706 全排列问题
// 按照字典序输出自然数 1 到 n 所有不重复的排列，即 n 的全排列，要求所产生的任一数字序列中不允许出现重复的数字。
#include <bits/stdc++.h>

inline bool isprime(int x) { 
	if (x == 1) return false; // 注意这步特判是必需的
	for (int i = 2; i * i <= x; ++i)
		if (x % i == 0)
			return false;
	return true;
}

const int N = 25;
int a[N], ans, n, k;

void dfs(int now, int sum, int sid) {
	// 现在已经选了 now 个数，当前总和为 sum
	// sid 是这次选数的起始下标，即我们从 a[sid] 开始选数枚举
	if (now == k) {
		if (isprime(sum))
			++ans;
		return ;
	}
	for (int i = sid; i <= n - k + now + 1; ++i)
		dfs(now + 1, sum + a[i], i + 1);
	return ;
}

int main() {
	scanf("%d%d", &n, &k);
	for (int i = 1; i <= n; ++i)
		scanf("%d", &a[i]);
	dfs(0, 0, 1);
	printf("%d\n", ans);
	return 0;
}










// 混合牛奶，由于乳制品产业利润很低，所以降低原材料（牛奶）价格就变得十分重要。帮助 Marry 乳业找到最优的牛奶采购方案。
#include<algorithm>
#include<iostream>
using namespace std;
struct farm {
	int price;
	int count;
}a[5005];
bool compare(const farm& a, const farm& b) {
	return a.price < b.price;
}
long long n, m;
int main()
{
	cin >> n >> m;
	long long res = 0;
	for (int i = 1; i <= m; i++) cin >> a[i].price>> a[i].count;
	sort(a + 1, a + 1 + m, compare);
	for (int i = 1; i <= m; i++) {
		if (n >= a[i].count) {
			res += a[i].price * a[i].count;
			n =n-a[i].count;
		}
		else {
			res += a[i].price * n;
			break;
		}
	}
	cout << res;
 
	return 0;
}

// 田忌赛马


#include <bits/stdc++.h>
using namespace std;
 
const int maxn=2005;
int t[maxn],q[maxn];
int money,tlow,qlow,thigh,qhigh;
 
int main() {
	int n;
	cin>>n;
	for(int i=0; i<n; i++) cin>>t[i];
	for(int i=0; i<n; i++) cin>>q[i];
	sort(t,t+n);
	sort(q,q+n);
 
	tlow=qlow=0;
	thigh=qhigh=n-1;
 
	while(tlow<=thigh) {
		if(t[thigh]>q[qhigh]) {
			money+=200;
			thigh--;
			qhigh--;
		} else if(t[thigh]<q[qhigh]) {
			money-=200;
			tlow++;
			qhigh--;
		} else {
			if(t[tlow]>q[qlow]) {
				money+=200;
				tlow++;
				qlow++;
			} else {
				if(t[tlow]<q[qhigh]) money-=200;
				tlow++;
				qhigh--;
			}
		}
	}
	cout<<money<<endl;
 
	return 0;
}


// 数列分段 对于给定的一个长度为 N 的正整数数列 A  i ​  ，现要将其分成连续的若干段，并且每段和不超过 M（可以等于 M），问最少能将其分成多少段使得满足要求。


#include<iostream>
using namespace std;
int main()
{
    int n,m,a[100000];
    int temp=0,result=1;
    int i;
    
    cin>>n>>m;
    for(i=1;i<=n;i++)	cin>>a[i];
    
    for(i=1;i<=n;i++)
    {
        if(a[i]+temp>m)	
        {
            result++;
            temp=0;
        } 
        temp+=a[i];
    }
    
    cout<<result<<endl;
    return 0;
}

// 纪念品分组 元旦快到了，校学生会让乐乐负责新年晚会的纪念品发放工作。为使得参加晚会的同学所获得的纪念品价值相对均衡，他要把购来的纪念品根据价格进行分组，但每组最多只能包括两件纪念品，并且每组纪念品的价格之和不能超过一个给定的整数。为了保证在尽量短的时间内发完所有纪念品，乐乐希望分组的数目最少。


#include <iostream>
#include <vector>
#include <algorithm>
 
int main() {
    int w, n;
    std::cin >> w;
    std::cin >> n;
 
    std::vector<int> prices(n);
    for (int i = 0; i < n; ++i) {
        std::cin >> prices[i];
    }
 
    std::sort(prices.begin(), prices.end());
 
    int left = 0, right = n - 1;
    int groups = 0;
 
    while (left <= right) {
        if (prices[left] + prices[right] <= w) {
            ++left;
            --right;
        } else {
            --right;
        }
        ++groups;
    } 
    std::cout << groups << std::endl;
 
    return 0;
}


// 合并果子，在一个果园里，多多已经将所有的果子打了下来，而且按果子的不同种类分成了不同的堆。多多决定把所有的果子合成一堆。

#include<iostream>
#include<algorithm>
#include<queue>
using namespace std;

long long sum = 0;  
priority_queue<int, vector<int>, greater<int>> q;  

bool comp(int a, int b) {
    return a < b;
}

int main() {
    int n;
    cin >> n;  

    for (int i = 1; i <= n; i++) {
        int x;
        cin >> x;
        q.push(x);  
    }

    while (q.size() != 1) {
        int a = q.top();  
        q.pop();
        int b = q.top(); 
        q.pop();

        sum += (a + b);  
        q.push(a + b);   
    }

    cout << sum << endl;  
    return 0;
}







// 装箱子，有一个箱子容量为 V，同时有 n 个物品，每个物品有一个体积。  现在从 n 个物品中，任取若干个装入箱内（也可以不取），使箱子的剩余空间最小。输出这个最小值。

#include<iostream>
using namespace std;
using ll = long long;
const int N = 2e4 + 9;
ll dp[N], v[N];
 
void solve()
{
	ll V, n; cin >> V >> n;
	for (int i = 1; i <= n; ++i) cin >> v[i];
	for (int i = 1; i <= n; ++i)
	{
		for (int j = V; j >= v[i]; --j)
		{
			dp[j] = max(dp[j], dp[j - v[i]] + v[i]);
		}
	}
 
	cout << V - dp[V];
}
 
int main()
{
	ios::sync_with_stdio(0), cin.tie(0), cout.tie(0);
	solve();
	return 0;
}



// 采药， 辰辰是个天资聪颖的孩子，他的梦想是成为世界上最伟大的医师。为此，他想拜附近最有威望的医师为师。医师为了判断他的资质，给他出了一个难题。医师把他带到一个到处都是草药的山洞里对他说：“孩子，这个山洞里有一些不同的草药

#include <iostream>

using namespace std;

const int N = 110, M = 1010;

int dp[M], t[N], w[N];
int n, m;

void input() {
    cin >> m >> n;
    for(int i = 1; i <= n; ++ i) cin >> t[i] >> w[i];
}

void solve() {
    for(int  i = 1; i <= n; ++ i) 
        for(int j = m; j >= t[i]; -- j) 
            dp[j] = max(dp[j], dp[j - t[i]] + w[i]);
}

int main(){
    input();

    solve();
    cout << dp[m] << endl;
    return 0;
}


// 开心金明，金明今天很开心，家里购置的新房就要领钥匙了，新房里有一间他自己专用的很宽敞的房间。更让他高兴的是，妈妈昨天对他说：“你的房间需要购买哪些物品

#include <bits/stdc++.h>
using namespace std;
const int M = 3e5+8, N = 3e3+8;
int f[M],w[N],c[N],m,n;
 
int main(void) {
	cin>>m>>n;
	for(int i=1; i <= n; i++){
		cin>>w[i]>>c[i];
		c[i] *= w[i];
	}
	for(int i=1; i <= n; i++) {
		for(int j=m; j >= w[i]; j--) {
			f[j] = max(f[j],c[i] + f[j - w[i]]);
		}
	}
	cout<<f[m];
	return 0;
}

// 能量项链，在 Mars 星球上，每个 Mars 人都随身佩带着一串能量项链。在项链上有 N 颗能量珠。能量珠是一颗有头标记与尾标记的珠子，这些标记对应着某个正整数。并且，对于相邻的两颗珠子，前一颗珠子的尾标记一定等于后一颗珠子的头标记。因为只有这样，

#include<bits/stdc++.h>
using namespace std;
int head[205],tail[205],f[205][205]={0};
int main()
{
	int ans=0,n,i,t,j,k;
	cin>>n;
	for(i=1;i<=n;i++) 
	{
		cin>>head[i];
		head[i+n]=head[i];
	}
	for(i=1;i<=2*n-1;i++) tail[i]=head[i+1];
	tail[2*n]=head[1];
	for(i=1;i<=2*n-1;i++) f[i][i]=0;
	for(t=1;t<=n-1;t++)
	{
		for(i=1;i<=2*n-t;i++)
		{
			j=i+t;
			for(k=i;k<=j-1;k++)
			{
				f[i][j]=max(f[i][j],f[i][k]+f[k+1][j]+head[i]*tail[k]*tail[j]);
			}
		}
	}
	for(i=1;i<=n;i++) ans=max(ans,f[i][i+n-1]);
	cout<<ans<<endl;
	return 0;
 }


// B2064 斐波那契数列

#include<stdio.h>
int main()
{
	int n,i,m,j;
	int x[40];
	x[1]=1;x[2]=1;
	scanf("%d",&n);
	for(i=1;i<=n;i++)
	{
		scanf("%d",&m);
		for(j=3;j<=m;j++)
		{
			x[j]=x[j-2]+x[j-1];
		}
		printf("%d\n",x[m]);
	}
}

//  一元三次方程求解

#include <bits/stdc++.h>
using namespace std;
double a,b,c,d;
int main()
{
    cin>>a>>b>>c>>d;
    for(double i=-100;i<=100;i+=0.001)
    {
        double x=i,y=x+0.001;
        double xx=x*x*x*a+x*x*b+x*c+d,yy=y*y*y*a+y*y*b+y*c+d;
        if(xx<=0&&yy>=0||xx>=0&&yy<=0)
        {
            double ans=(x+y)/2;
            cout<<fixed<<setprecision(2)<<ans<<" ";
        }
    }
      return 0;
}

// 字符序列的子序列是指从给定字符序列中随意地（不一定连续）去掉若干个字符（可能一个也不去掉）后所形成的字符序列。令给定的字符序列 X={x  0 ​  ,x  1 ​  ,⋯,x  m−1 ​  }，

// Code by KSkun, 2018/6
#include <cstdio>
#include <cctype>
#include <cstring>

#include <algorithm>

const int MAXN = 5005, MO = 1e8;

int n, m, dp[2][MAXN], dp2[2][MAXN];
char s1[MAXN], s2[MAXN];

int main() {
    scanf("%s%s", s1 + 1, s2 + 1);
    n = strlen(s1 + 1) - 1; m = strlen(s2 + 1) - 1;
    for(int i = 0; i <= m; i++) dp2[0][i] = 1;
    dp2[1][0] = 1;
    int p = 1;
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= m; j++) {
            dp[p][j] = std::max(dp[p ^ 1][j], dp[p][j - 1]);
            dp2[p][j] = 0;
            if(s1[i] == s2[j]) {
                dp[p][j] = std::max(dp[p][j], dp[p ^ 1][j - 1] + 1);
            }
            if(s1[i] == s2[j] && dp[p][j] == dp[p ^ 1][j - 1] + 1) {
                dp2[p][j] = (dp2[p][j] + dp2[p ^ 1][j - 1]) % MO;
            }
            if(dp[p][j] == dp[p ^ 1][j]) {
                dp2[p][j] = (dp2[p][j] + dp2[p ^ 1][j]) % MO;
            }
            if(dp[p][j] == dp[p][j - 1]) {
                dp2[p][j] = (dp2[p][j] + dp2[p][j - 1]) % MO;
            }
            if(dp[p][j] == dp[p ^ 1][j - 1]) {
                dp2[p][j] = ((dp2[p][j] - dp2[p ^ 1][j - 1]) % MO + MO) % MO;
            }
        }
        p ^= 1;
    }
    printf("%d\n%d", dp[p ^ 1][m], dp2[p ^ 1][m]);
    return 0;
}

// 在一个圆形操场的四周摆放 N 堆石子，现要将石子有次序地合并成一堆，规定每次只能选相邻的 2 堆合并成新的一堆，并将新的一堆的石子数，记为该次合并的得分。

#include<iostream>
#include<cstring>
using namespace std;
const int N=500+5;
int s[N];
int Max[N][N];
int Min[N][N];
int MAX=-1;
int MIN=1e9;
int main(){
    int n;
    cin>>n;
    memset(Max,-0x3f,sizeof(Max));
    memset(Min,0x3f,sizeof(Min));
    for(int i=1;i<=n;i++){
        cin>>s[i];
        s[i+n]=s[i];
        Max[i][i]=0;
        Max[i+n][i+n]=0;
        Min[i][i]=0;
        Min[i+n][i+n]=0;
    }
    for(int i=1;i<=2*n;i++)s[i]+=s[i-1];
    for(int i=1;i<=n;i++){
        for(int len=2;len<=n;len++){
            for(int l=i;l+len-1<=n+i-1;l++){
                int r=l+len-1;
                for(int k=l;k<r;k++){
                    Max[l][r]=max(Max[l][r],Max[l][k]+Max[k+1][r]+s[r]-s[l-1]);
                    Min[l][r]=min(Min[l][r],Min[l][k]+Min[k+1][r]+s[r]-s[l-1]);
                }
               
            }
        }
        MAX=max(MAX,Max[i][n+i-1]);
        MIN=min(MIN,Min[i][n+i-1]);
    }
    cout<<MIN<<endl<<MAX;
    return 0;
}


// 给定非负整数 n，求 2  n   的值。


#include <bits/stdc++.h>
using namespace std;
int main() {
	long long n;
	cin>>n;
	int t=1;
	for(int i=1;i<=n;i++){
		t*=2;
	}
	cout<<t; 
	return 0;
}
```
