from pydantic import BaseModel


class GenomePosition(BaseModel):
    chrom: str
    pos: int

    def __hash__(self) -> int:
        return (self.chrom,self.pos).__hash__()
    
    def __lt__(self, other) -> bool:
        if self.chrom != other.chrom:
            raise ValueError(f"Cannot compare positions on different chromosomes: {self.chrom} and {other.chrom}")
        return self.pos < other.pos

class GenomeRange(BaseModel):
    chrom: str
    start: int
    end: int

    def __contains__(self,pos: GenomePosition):
        if pos.chrom==self.chrom and pos.pos >= self.start and pos.pos <= self.end:
            return True
        else:
            return False
