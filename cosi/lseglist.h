
namespace cosi {
namespace lseglist {

typedef double loc_t;

class LazyLoc {
public:
	 LazyLoc(): locComputed( false ) { }
	 loc_t getLoc() const { if ( !locComputed ) loc = computeLoc(); return loc; }
protected:
	 loc_t computeLoc() const = 0;
private:
	 mutable loc_t loc;
	 mutable bool locComputed;
};
typedef boost::shared_ptr<LazyLoc> LazyLocP;

class LazySeglist {

	 LazyLocP beg;
	 LazyLocP end;
	 
	 virtual loc_t getBeg() const = 0;
	 virtual loc_t getEnd() const = 0;

	 virtual std::pair<LSeglistP,LSeglistP> split( loc_t loc ) = 0;
};

class LSeglist_OneSeg: public LSeglist {

	 LSeglist_OneSeg( loc_t beg_, loc_t end_ ): beg( beg_ ), end( end_ ) { }
	 
	 virtual loc_t getBeg() const { return beg; }
	 virtual loc_t getEnd() const { return end; }

	 virtual std::pair<LSeglistP,LSeglistP> split( loc_t loc ) {
		 LSeglistP right = boost::make_shared<LSeglist_OneSeg( loc, end );
		 this->end = loc;
		 return std::make_pair( shared_from_this(), right );
	 }
	 
private:
	 loc_t beg, end;
};

class LSeglist_CachedSeg: public LSeglist {

public:	 
	 virtual loc_t getBeg() const { if ( !begComputed ) beg = computeBeg(); return beg; }
	 virtual loc_t getEnd() const { if ( !endComputed ) end = computeEnd(); return end; }
	 
protected:
	 LSeglist_CachedSeg(): begComputed( false ), endComputed( false ) { }

	 virtual loc_t computeBeg() const = 0;
	 virtual loc_t computeEnd() const = 0;
	 
private:
	 mutable loc_t beg, end;
	 mutable bool_t begComputed, endComputed;
};

class LSeglist_Union: public LSeglist_CachedSeg {

protected:
	 virtual loc_t computeBeg() const { return boost::min( lseglist1->getBeg(), lseglist2->getBeg() ); }
	 virtual loc_t computeEnd() const { return boost::max( lseglist1->getEnd(), lseglist2->getEnd() ); }

	 virtual std::pair<LSeglistP,LSeglistP> split( loc_t loc );
	 
private:
	 LSeglistP lseglist1, lseglist2;
	 
};

class LSeglist_Union_SplitLeft: public LSeglist_CachedSeg {

protected:
	 virtual loc_t computeBeg() const { return lseglistUnion->lseglist1->computeBeg(); }
	 virtual loc_t computeEnd() const {
		 
	 }
	 
private:
	 LSeglist_UnionP lseglistUnion;
};

}
}
