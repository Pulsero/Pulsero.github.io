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


```
