// [[Rcpp::depends(IRanges,S4Vectors)]]

#include <Rcpp.h>
#include <Rdefines.h>

#include <IRanges_interface.h>
#include <_IRanges_stubs.c>
// #include <tuple> // need support C++11
// #include <S4Vectors_interface.h>

// #include "IntervalTree.h"

using namespace std;
using namespace Rcpp;


// IntervalTree.h

#ifndef __INTERVAL_TREE_H
#define __INTERVAL_TREE_H

#include <vector>
#include <algorithm>
#include <iostream>


template <class T, typename K = int>
class Interval {
public:
    K start;
    K stop;
    T value;
    Interval(K s, K e, const T& v)
        : start(s)
        , stop(e)
        , value(v)
    { }
};


class IntervalRange {
public:
  int start;
  int stop;
    IntervalRange(int s, int e)
      : start(s)
      , stop(e)
  { }
};



template <class T, typename K>
int intervalStart(const Interval<T,K>& i) {
    return i.start;
}

template <class T, typename K>
int intervalStop(const Interval<T,K>& i) {
    return i.stop;
}

template <class T, typename K>
ostream& operator<<(ostream& out, Interval<T,K>& i) {
    out << "Interval(" << i.start << ", " << i.stop << "): " << i.value;
    return out;
}

template <class T, typename K = int>
class IntervalStartSorter {
public:
    bool operator() (const Interval<T,K>& a, const Interval<T,K>& b) {
        return a.start < b.start;
    }
};


// IntervalTree class

template <class T, typename K = int>
class IntervalTree {
  
public:
  typedef Interval<T,K> interval;
  typedef vector<interval> intervalVector;
  typedef IntervalTree<T,K> intervalTree;
  
  intervalVector intervals;
  intervalTree* left;
  intervalTree* right;
  int center;
  
  IntervalTree<T,K>(void)
  : left(NULL)
  , right(NULL)
  , center(0)
  { }
  
  IntervalTree<T,K>(const intervalTree& other) {
    center = other.center;
    intervals = other.intervals;
    if (other.left) {
      left = (intervalTree*) malloc(sizeof(intervalTree));
      *left = *other.left;
    } else {
      left = NULL;
    }
    if (other.right) {
      right = new intervalTree();
      *right = *other.right;
    } else {
      right = NULL;
    }
  }
  
  IntervalTree<T,K>& operator=(const intervalTree& other) {
    center = other.center;
    intervals = other.intervals;
    if (other.left) {
      left = new intervalTree();
      *left = *other.left;
    } else {
      left = NULL;
    }
    if (other.right) {
      right = new intervalTree();
      *right = *other.right;
    } else {
      right = NULL;
    }
    return *this;
  }
  
  IntervalTree<T,K>(
		    intervalVector& ivals,
		    unsigned int depth = 16,
		    unsigned int minbucket = 64,
		    int leftextent = 0,
		    int rightextent = 0,
		    unsigned int maxbucket = 512
		    )
  : left(NULL)
    , right(NULL)
  {
    
    --depth;
    IntervalStartSorter<T,K> intervalStartSorter;
    if (depth == 0 || (ivals.size() < minbucket && ivals.size() < maxbucket)) {
      sort(ivals.begin(), ivals.end(), intervalStartSorter);
      intervals = ivals;
    } else {
      if (leftextent == 0 && rightextent == 0) {
	// sort intervals by start
	sort(ivals.begin(), ivals.end(), intervalStartSorter);
      }
      
      int leftp = 0;
      int rightp = 0;
      int centerp = 0;
      
      if (leftextent || rightextent) {
	leftp = leftextent;
	rightp = rightextent;
      } else {
	leftp = ivals.front().start;
	vector<K> stops;
	stops.resize(ivals.size());
	transform(ivals.begin(), ivals.end(), stops.begin(), intervalStop<T,K>);
	rightp = *max_element(stops.begin(), stops.end());
      }
      
      //centerp = ( leftp + rightp ) / 2;
      centerp = ivals.at(ivals.size() / 2).start;
      center = centerp;
      
      intervalVector lefts;
      intervalVector rights;
      
      for (typename intervalVector::iterator i = ivals.begin(); i != ivals.end(); ++i) {
	interval& interval = *i;
	if (interval.stop < center) {
	  lefts.push_back(interval);
	} else if (interval.start > center) {
	  rights.push_back(interval);
	} else {
	  intervals.push_back(interval);
	}
      }
      
      if (!lefts.empty()) {
	left = new intervalTree(lefts, depth, minbucket, leftp, centerp);
      }
      if (!rights.empty()) {
	right = new intervalTree(rights, depth, minbucket, centerp, rightp);
      }
    }
  }
  
  void findWithin(K start, K stop, map<string,string>& contained) {
    if (!intervals.empty() && ! (stop < intervals.front().start)) {
      for (typename intervalVector::iterator i = intervals.begin(); i != intervals.end(); ++i) {
	interval& interval = *i;
	if (start >= interval.start && stop <= interval.stop) {
	  vector<string> insVal = interval.value;
	  contained.insert(map<string,string>::value_type(insVal[0],insVal[1]));
	}
      }
    }
    
    if (left && start <= center) {
      left->findWithin(start, stop, contained);
    }
    
    if (right && stop >= center) {
      right->findWithin(start, stop, contained);
    }
    
  }
  
  // find query to have same start as tree interval
  void findStart(K start, K stop, map<string,string>& contained) {
    if (!intervals.empty() && ! (stop < intervals.front().start)) {
      for (typename intervalVector::iterator i = intervals.begin(); i != intervals.end(); ++i) {
	interval& interval = *i;
	if (start == interval.start && stop <= interval.stop) {
	  vector<string> insVal = interval.value;
	  contained.insert(map<string,string>::value_type(insVal[0],insVal[1]));
	}
      }
    }
    
    if (left && start <= center) {
      left->findStart(start, stop, contained);
    }
    
    if (right && stop >= center) {
      right->findStart(start, stop, contained);
    }
  }
  
  

  // find query to have same end as tree interval
  void findEnd(K start, K stop, map<string,string>& contained) {
    if (!intervals.empty() && ! (stop < intervals.front().start)) {
      for (typename intervalVector::iterator i = intervals.begin(); i != intervals.end(); ++i) {
	interval& interval = *i;
	if (start >= interval.start && stop == interval.stop) {
	  vector<string> insVal = interval.value;
	  contained.insert(map<string,string>::value_type(insVal[0],insVal[1]));
	}
      }
    }
    
    if (left && start <= center) {
      left->findEnd(start, stop, contained);
    }
    
    if (right && stop >= center) {
      right->findEnd(start, stop, contained);
    }
  }

  // find query to have exact overlap as tree interval
  void findExact(K start, K stop, map<string,string>& contained) {
    if (!intervals.empty() && ! (stop < intervals.front().start)) {
      for (typename intervalVector::iterator i = intervals.begin(); i != intervals.end(); ++i) {
	interval& interval = *i;
	if (start == interval.start && stop == interval.stop) {
	  vector<string> insVal = interval.value;
	  contained.insert(map<string,string>::value_type(insVal[0],insVal[1]));
	}
      }
    }
    
    if (left && start <= center) {
      left->findExact(start, stop, contained);
    }
    
    if (right && stop >= center) {
      right->findExact(start, stop, contained);
    }
  }
  
  ~IntervalTree(void) {
    // traverse the left and right
    // delete them all the way down
    if (left) {
      delete left;
    }
    if (right) {
      delete right;
    }
  }
  
};



// IntervalForest Class: should hold as many trees as chromosomes. Use a map to retrive them with 'string'

// template <class T, typename K = int>
// class IntervalForest
// {
// public:
//   typedef IntervalTree<T,K> intervalTree;
//   map<string,intervalTree> forest;

// }



#endif

///////////////////////////////////////////////////////////////

typedef union {void *p; DL_FUNC fn;} fn_ptr;

DL_FUNC R_ExternalPtrAddrFn(SEXP s){
     fn_ptr tmp;
     tmp.p =  EXTPTR_PTR(s);
     return tmp.fn;
}



typedef Interval< vector<string> > interval;
typedef Interval<string> iRange;
typedef vector<interval> intervalVector;
typedef vector<iRange> rangeVector;
typedef IntervalTree< vector<string> > intervalTree;
typedef IntegerVector::iterator intIter;
typedef CharacterVector::iterator charIter;


/*
  makeTree

  INPUT: unlistData slot from a GRangeList
  OUTPUT: pointer to a IntervalTree object

 */


// [[Rcpp::export]]
SEXP makeTree(SEXP r_unlistData)
{
  //  BEGIN_RCPP
  intervalVector intervals;

  // We need to partition within seqnames!!!
  // for the test we just use chr4
  S4 unlistData(r_unlistData);

  // S4 seqRle = unlistData.slot("seqnames");
  // IntegerVector seqVals = seqRle.slot("values");
  // vector< string > seqLev = as< vector< string > >(seqVals.attr("levels"));
  
 
  S4 strandRle = unlistData.slot("strand");
  IntegerVector strandVal = strandRle.slot("values");
  IntegerVector strandLen = strandRle.slot("lengths");

  vector< int > strandValues;

  for(pair<intIter, intIter> itr(strandLen.begin(),strandVal.begin()); itr.first != strandLen.end(); 
      ++itr.first, ++itr.second)
    for(int i = 0; i != *itr.first; ++i )
      {
      strandValues.push_back(*itr.second);
    }

  vector< string > strandLevs = as< vector< string > >(strandVal.attr("levels"));

  DataFrame metadata = unlistData.slot("elementMetadata");
  S4 ranges = unlistData.slot("ranges");
  
  IntegerVector start = ranges.slot("start");
  IntegerVector width = ranges.slot("width");
  vector<string> exname = as< vector<string>  >(metadata["exon_name"]);
  vector<string> txname = as< vector<string>  >(metadata["tx_name"]);

  for(size_t i; i != start.size();++i)
    {
      int stop = start(i) + width(i) - 1;
      int id = strandValues[i] - 1;
      vector<string> vec{txname[i], exname[i], strandLevs[id]};
      intervals.push_back(interval(start(i),stop,vec));
    }


  Rcpp::XPtr< intervalTree > tree_ptr(new intervalTree,true);
  *tree_ptr =  intervalTree(intervals);

  return(tree_ptr);
  // END_RCPP
}



// Templated function to map match results from right & left mates.
// Key = tx_name
// Values = ex_names

template<typename KeyType, typename LeftValue, typename RightValue>
map<KeyType, pair<LeftValue, RightValue> > IntersectMaps(const map<KeyType, LeftValue> & left, 
							 const map<KeyType, RightValue> & right)
{
  map<KeyType, pair<LeftValue, RightValue> > result;
  typename map<KeyType, LeftValue>::const_iterator il = left.begin();
  typename map<KeyType, RightValue>::const_iterator ir = right.begin();
  while (il != left.end() && ir != right.end())
    {
      if (il->first < ir->first)
	++il;
      else if (ir->first < il->first)
	++ir;
      else
        {
	  result.insert(make_pair(il->first, make_pair(il->second, ir->second)));
	  ++il;
	  ++ir;
        }
    }
  return result;
}


// find overlaps for gapped reads
map<string,string> matchGapped(IRanges_holder* irange, string strand, intervalTree& forest, int lenElt)
{

  int id = 0;
  map<string,string> results;
  int i_start, i_width, i_end;
  rangeVector queries;
  
  // create a vector of intervals. They all have to map to the same 'junction' exon
  // junction exons are coded as 2 or more intervals, same ex_name
  for(int id=0; id != lenElt; ++id)
    {
      i_start = get_start_elt_from_IRanges_holder(irange,id);
      i_width = get_width_elt_from_IRanges_holder(irange,id);
      i_end = i_start + i_width - 1;
      
      
      queries.push_back(iRange(i_start,i_end,strand));
    }


  // We expect ALL element of a gapped range to overlap EXACTLY to the same tx-exon pairs
  // We check this iteratively setting the first (the left-most element on + strand) overlap as pivot:
  // all others must be equal to that

  for (rangeVector::iterator q = queries.begin(); q != queries.end(); ++q) {
    map<string,string> cache_search;
    if(q == queries.begin())
      {
	forest.findEnd(q->start, q->stop, cache_search);

	if(cache_search.size() == 0) // no match...stop!
	  return(results);
	else
	  results = cache_search;
      }
    else if(q < (queries.end()-1)) // all the chunks but first and last
      {
	forest.findWithin(q->start, q->stop, cache_search);
	
	if(cache_search.size() == 0 || cache_search != results)
	  return(map<string,string>()); // no match here
      }
    else {

      forest.findStart(q->start, q->stop, cache_search);

      if(cache_search.size() == 0 || cache_search != results)
	return(map<string,string>()); // no match here
    }
  }

  // all 'chunks' must map to the same exon (in Sequgio's txdb junction regions are coded as a single region)
  return(results);
  
}

// find overlaps for ungapped reads

map<string,string> matchSimple(IRanges_holder* irange, string strand, intervalTree& forest)
{
  int id = 0;
  int i_start = get_start_elt_from_IRanges_holder(irange,id);
  int i_width = get_width_elt_from_IRanges_holder(irange,id);
  int i_end = i_start + i_width - 1;
  
  rangeVector queries;
  queries.push_back(iRange(i_start,i_end,strand));

  map<string,string> results;
  
  for (rangeVector::iterator q = queries.begin(); q != queries.end(); ++q) {
    forest.findWithin(q->start, q->stop, results);
  }

  return(results);

}


// [[Rcpp::export]]
SEXP getOverlaps(SEXP r_forest, SEXP r_reads)
{
  BEGIN_RCPP


  // cigar_ranges ///////////////////////////////////////////////////////////////////////////

  SEXP xp;
  List nativeSymbolInfo;
  Function getNativeSymbolInfo("getNativeSymbolInfo");
  
  nativeSymbolInfo = getNativeSymbolInfo("cigar_ranges");
  xp = nativeSymbolInfo["address"];
  DL_FUNC  cigar_ranges_p = R_ExternalPtrAddrFn( xp ) ;
  static SEXP (*cigar_ranges)(SEXP cigar, SEXP flag, SEXP space, SEXP pos, SEXP f,
  			      SEXP ops, SEXP drop_empty_ranges, SEXP reduce_ranges,
  			      SEXP with_ops) = NULL;
  cigar_ranges = (SEXP(*)(SEXP, SEXP,SEXP, SEXP,SEXP, SEXP,SEXP, SEXP,SEXP)) cigar_ranges_p;

  ///////////////////////////////////////////////////////////////////////////////////////////

  S4 forest_obj(r_forest);
  SEXP forest_ptr;
  forest_ptr = forest_obj.slot("ptr");
  
  intervalTree *forest_p = static_cast<intervalTree*>(R_ExternalPtrAddr(forest_ptr));
  intervalTree&  forest = *forest_p;

  List reads(r_reads);

  IntegerVector groupid = reads["groupid"];
  IntegerVector strand_v = reads["strand"]; //1=+; 2=-; *=3
  CharacterVector cigar = reads["cigar"];
  int nreads = groupid.size();


  const vector<string> strandLev{"+","-","*"}; 
  vector < string > strand;
  
  for(intIter itr = strand_v.begin(); itr != strand_v.end(); ++itr)
    {
      string tmpVal = strandLev[(*itr)-1];
      strand.push_back(tmpVal);
    }
  

  
  Environment GA("package:GenomicAlignments");
  CharacterVector CIGAR_OPS = GA["CIGAR_OPS"];

  SEXP flag = R_NilValue;
  SEXP f = R_NilValue;
  LogicalVector with_ops(1), drop_empty_ranges(1), reduce_ranges(1);
  with_ops[0] = 0;
  drop_empty_ranges[0] = 1;
  reduce_ranges[0] = 0;

  IntegerVector pos, space(1);
  pos = reads["pos"];

  space[0] = 3;


  // "CompressedIRangesList": every element is a IRanges object
  SEXP cigar_rng = (*cigar_ranges)(cigar, flag, space, pos, f, CIGAR_OPS, drop_empty_ranges, 
  				reduce_ranges,with_ops);

  CompressedIRangesList_holder cigar_hold = hold_CompressedIRangesList(cigar_rng);


  vector<int> lenOut;

  // Loop along reads //


  int i_start, i_width, i_end;

  // unordered set or set?
  vector< string > collect_results;

  for(size_t i=0; i != nreads; i+=2) // increment by 2 as reads are paired
    {

      string i_strand = strand[i];
      map<string,string> result_plus, result_neg;
      int len_elt;
      IRanges_holder elt;

      // exon names pairs in counting are listed in direction + -> -, so we always start from "+" mate
      int j = (i_strand == "+") ? i : (i+1); // recall: condition ? value_if_true : value_if_false

      // This should go in a function that process the i-th read

      // place 1st mate ("+")
      elt = get_elt_from_CompressedIRangesList_holder(&cigar_hold,j);
      len_elt = get_length_from_IRanges_holder(&elt);
      
      if(len_elt == 1)
	{
	  result_plus = matchSimple(&elt, i_strand,forest);
	  // some reads might fall in exons that are shared by different tx, so appear more than once!
	  // if(result_plus.size() > 1)
	  //   {
	  //     copy(result_plus.begin(), result_plus.end(), ostream_iterator<string>(cout, "-"));
	  //     cout << j << endl;
	  //   }
	}
      else
	{
	  // matchGapped
	}
 
      
      // if no overlap for + read
      if(result_plus.size() == 0)
      	continue;

      // get next read, - strand. Is it listed before or after + strand? -> k
      int k = (j > i) ? i : i+1;

      elt = get_elt_from_CompressedIRangesList_holder(&cigar_hold,k);
      len_elt = get_length_from_IRanges_holder(&elt);
      
      if(len_elt == 1)
	{
	  result_neg = matchSimple(&elt, strand[k],forest);
	}
      else
	{
	  result_neg = matchGapped(&elt, strand[k],forest,len_elt);
	}

      
      // no match from - strand read
      if(result_neg.size() == 0)
      	continue;

      // cout << i << " | " << result_plus.size()<< "/" << result_neg.size();

      map<string, pair<string,string> > merge_results = IntersectMaps(result_plus,result_neg);
      // we might get no intersection!!!
      // cout << "=> " << merge_results.size() << endl;
      
      if(merge_results.size() > 0)
	{
	  for(map<string, pair<string,string> >::iterator iter=merge_results.begin();
	      iter != merge_results.end(); ++iter) {

	    string ss = iter->second.first + "." + iter->second.second;
	    collect_results.push_back(ss);
	  }
	  
	}

    }


  int ans = 0;
  return(wrap(collect_results));
  
  END_RCPP
}

